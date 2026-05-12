"""Regenerate numeric reference snapshots for tests/test_numeric_snapshots.py.

For each example model and each solver method, this script regenerates two
kinds of reference data and stores them as ``.npz`` files in
``tests/numeric_refs/``:

Static fingerprint (one ``<model>_static.npz`` per model)
    ``q0``, ``q0_dot``  -- initial coords / velocities packed by the solver.
    ``M``               -- diagonal mass matrix at t = 0.
    ``Phi_q``           -- constraint Jacobian at t = 0.
    ``gamma``           -- acceleration-level RHS at t = 0.
    ``Q``               -- generalized-force vector at t = 0.

Dynamic fingerprint (one ``<model>_<solver>.npz`` per (model, solver) pair)
    ``t_samples``       -- 11 equally-spaced time stamps (decile indices).
    ``u_samples``       -- corresponding rows of ``uT``.

CLI
---
    python tests/_regen_refs.py                       # all models, both solvers
    python tests/_regen_refs.py --model msd           # one model, both solvers
    python tests/_regen_refs.py --solver lsoda        # all models, one solver
    python tests/_regen_refs.py --model cs --solver casadi
"""
from __future__ import annotations

import argparse
import importlib
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
PMD_ROOT = HERE.parent
REFS_DIR = HERE / "numeric_refs"
REFS_DIR.mkdir(exist_ok=True)

# Make `examples` importable when running this script directly.
if str(PMD_ROOT) not in sys.path:
    sys.path.insert(0, str(PMD_ROOT))

# ---------------------------------------------------------------------------
# Model registry
# ---------------------------------------------------------------------------
MODELS: dict[str, str] = {
    "msd": "examples.mass_spring_damper",
    "cs":  "examples.crank_slider",
    "fbl": "examples.four_bar_linkage",
}

SOLVERS: tuple[str, ...] = ("lsoda", "casadi")

# Solver method strings expected by ``PlanarMultibodyModel.solve``.
_METHOD = {
    "lsoda":  "LSODA",
    "casadi": "CASADI-DAE",
}

# Number of dynamic samples (decile indices: 0%, 10%, ..., 100%).
N_SAMPLES = 11


# ---------------------------------------------------------------------------
# Fingerprint helpers
# ---------------------------------------------------------------------------
def static_fingerprint(model) -> dict[str, np.ndarray]:
    """Snapshot the constant matrices and the t=0 evaluations of Phi_q, gamma, Q."""
    # `_initialize` already ran in PlanarMultibodyModel.__init__, so the
    # marker fields are populated for the initial pose. We can call the
    # private hooks directly. Time-dependent driver functions need ``model.t``.
    model.t = 0.0
    q0 = model._bodies2u()
    Phi_q = model._compute_jacobian()
    gamma = model._rhs_acceleration()
    Q = model._compute_force()
    return {
        "q0":     np.ascontiguousarray(q0),
        "M":      np.ascontiguousarray(model.M_matrix),
        "Phi_q":  np.ascontiguousarray(Phi_q),
        "gamma":  np.ascontiguousarray(gamma),
        "Q":      np.ascontiguousarray(Q),
    }


def dynamic_fingerprint(model, method: str, t_final: float, n_eval: int,
                        ic_correct: bool) -> dict[str, np.ndarray]:
    """Run a full simulation and sample 11 decile points from the trajectory."""
    T, uT = model.solve(
        method=method,
        t_final=t_final,
        t_eval=np.linspace(0.0, t_final, n_eval),
        ic_correct=ic_correct,
    )
    idx = np.linspace(0, len(T) - 1, N_SAMPLES).astype(int)
    return {
        "t_samples": np.ascontiguousarray(T[idx]),
        "u_samples": np.ascontiguousarray(uT[idx, :]),
    }


# ---------------------------------------------------------------------------
# Per-model driver
# ---------------------------------------------------------------------------
def regenerate(model_key: str, solvers: tuple[str, ...]) -> None:
    module_name = MODELS[model_key]
    print(f"\n=== {model_key} ({module_name}) ===")

    # Re-import fresh each call so that side-effecting Body() / Ground.add_marker
    # state from a previous regeneration cannot leak through.
    if module_name in sys.modules:
        del sys.modules[module_name]
    mod = importlib.import_module(module_name)

    # Static snapshot (build once, fingerprint, save)
    print("  static  ...", end=" ", flush=True)
    model = mod.build_model()
    static = static_fingerprint(model)
    static_path = REFS_DIR / f"{model_key}_static.npz"
    np.savez(static_path, **static)
    print(f"OK  ({static_path.name})")

    # Dynamic snapshot per solver (rebuild each time to avoid stale state)
    for solver in solvers:
        print(f"  {solver:<7} ...", end=" ", flush=True)
        model = mod.build_model()
        try:
            dyn = dynamic_fingerprint(
                model,
                method=_METHOD[solver],
                t_final=mod.T_FINAL,
                n_eval=mod.N_EVAL,
                ic_correct=mod.IC_CORRECT,
            )
        except Exception as exc:  # noqa: BLE001
            print(f"FAILED ({type(exc).__name__}: {exc})")
            continue
        dyn_path = REFS_DIR / f"{model_key}_{solver}.npz"
        np.savez(dyn_path, **dyn)
        print(f"OK  ({dyn_path.name})")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--model",
        choices=sorted(MODELS),
        default=None,
        help="Regenerate only this model (default: all).",
    )
    parser.add_argument(
        "--solver",
        choices=sorted(SOLVERS),
        default=None,
        help="Regenerate only this solver's dynamic refs (default: both).",
    )
    args = parser.parse_args()

    keys = [args.model] if args.model else list(MODELS)
    solvers = (args.solver,) if args.solver else SOLVERS
    for k in keys:
        regenerate(k, solvers)
    print("\nAll done.")


if __name__ == "__main__":
    main()
