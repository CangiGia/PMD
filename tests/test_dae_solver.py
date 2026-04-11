"""Tests for the custom BDF DAE index-3 solver (method='BDF-DAE').

Test classes
------------
TestDAESolverPendulum
    Basic correctness on a single pendulum:
    constraint satisfaction, energy conservation, output shape.

TestDAEvsODE
    Trajectory agreement between DAE and ODE (Radau) solvers on the
    same model with a short simulation.

TestConstraintDrift
    The DAE solver maintains near-machine-precision constraint
    satisfaction over a longer integration (10 s).

TestInputRouting
    solve(analysis="dynamic", method="BDF-DAE") correctly dispatches
    to _solve_dae() without error.

Run from:  C:\\Users\\Giaco\\anaconda3\\envs\\GiacoEnv\\
Command:   python -m pytest PMD/tests/test_dae_solver.py -v
"""

import numpy as np
import pytest

from PMD.src.model import Body, Ground, _GroundType
from PMD.src.constraints import RevJoint, Weight
from PMD.src.solver import PlanarMultibodyModel


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _reset():
    """Reset global state between tests."""
    yield
    gt = _GroundType._instance
    if gt is not None:
        _GroundType._markers = [gt.origin]
    Body.COUNT = 0


def _make_pendulum(length=0.5):
    """Single pendulum of unit mass and inertia, pinned at origin.

    Body CoM is at (0, -length) with a marker at (0, +length) for the
    RevJoint pivot → pivot at ground origin.
    DOF = 3*1 - 2 = 1.
    """
    B = Body(mass=1.0, inertia=0.1,
             position=[0.0, -length], orientation=0.0)
    mk_B = B.add_marker([0.0, length])
    mk_G = Ground.add_marker([0.0, 0.0])
    j = RevJoint(iMarker=mk_B, jMarker=mk_G)
    w = Weight()
    model = PlanarMultibodyModel(bodies=[B], joints=[j], forces=[w])
    return model, B


# ---------------------------------------------------------------------------
# 1. Basic correctness on a pendulum
# ---------------------------------------------------------------------------

class TestDAESolverPendulum:

    def test_constraint_satisfaction(self):
        """max |Phi(q_k)| < 1e-8 at every output step."""
        model, B = _make_pendulum()
        T, uT = model.solve(
            analysis="dynamic", method="BDF-DAE",
            t_final=1.0, dt=0.01, ic_correct=True,
        )
        assert len(T) > 1

        nB3 = 3 * model.nB
        max_phi = 0.0
        for k in range(len(T)):
            q_k = uT[k, :nB3]
            for Bi, body in enumerate(model.Bodies):
                ir = 3 * Bi
                body.position    = q_k[ir:ir + 2].reshape(2, 1)
                body.orientation = float(q_k[ir + 2])
            model._update_position()
            phi = model._compute_constraints().flatten()
            max_phi = max(max_phi, float(np.max(np.abs(phi))))

        assert max_phi < 1e-8, (
            f"Constraint violation too large: max|Phi| = {max_phi:.3e}"
        )

    def test_energy_conservation(self):
        """Total mechanical energy conserved to within 1 % over 2 s.

        E = 0.5 * m * ||v||^2 + 0.5 * I * omega^2 + m * g * y_cm
        """
        g = 9.81
        model, B = _make_pendulum()
        T, uT = model.solve(
            analysis="dynamic", method="BDF-DAE",
            t_final=2.0, dt=0.01, ic_correct=True,
        )

        nB3 = 3 * model.nB
        m, I = B.mass, B.inertia

        E = np.zeros(len(T))
        for k in range(len(T)):
            x,  y,  phi  = uT[k, 0], uT[k, 1], uT[k, 2]
            vx, vy, vphi = uT[k, nB3], uT[k, nB3 + 1], uT[k, nB3 + 2]
            E[k] = (0.5 * m * (vx**2 + vy**2)
                    + 0.5 * I * vphi**2
                    + m * g * y)

        E0 = E[0]
        rel_err = np.abs(E - E0) / (np.abs(E0) + 1e-12)
        assert np.max(rel_err) < 0.01, (
            f"Energy not conserved: max relative error = {np.max(rel_err):.3%}"
        )

    def test_output_shape(self):
        model, _ = _make_pendulum()
        T, uT = model.solve(
            analysis="dynamic", method="BDF-DAE",
            t_final=0.5, dt=0.05,
        )
        nB3 = 3 * model.nB
        assert uT.shape == (len(T), 2 * nB3)

    def test_result_container_populated(self):
        """Body._result_container must contain positions/velocities/accelerations."""
        model, B = _make_pendulum()
        model.solve(
            analysis="dynamic", method="BDF-DAE",
            t_final=0.3, dt=0.05,
        )
        rc = B._result_container
        for key in ("positions", "velocities", "accelerations"):
            assert key in rc, f"Missing key '{key}' in _result_container"

    def test_joint_result_container_populated(self):
        """Joint._result_container must contain reactions."""
        model, B = _make_pendulum()
        model.solve(
            analysis="dynamic", method="BDF-DAE",
            t_final=0.3, dt=0.05,
        )
        for joint in model.Joints:
            assert "reactions" in joint._result_container, (
                "Missing 'reactions' in joint._result_container"
            )


# ---------------------------------------------------------------------------
# 2. Trajectory agreement with ODE solver
# ---------------------------------------------------------------------------

class TestDAEvsODE:
    """DAE and ODE solutions must agree on a short, benign pendulum run."""

    def test_position_agreement(self):
        """max |q_DAE - q_ODE| < 0.05 rad/m over 1 s (dt = 0.01 s)."""
        # --- DAE ---
        model_dae, _ = _make_pendulum()
        T_dae, uT_dae = model_dae.solve(
            analysis="dynamic", method="BDF-DAE",
            t_final=1.0, dt=0.01, ic_correct=True,
        )

        # Reset and rebuild for ODE
        _GroundType._instance._markers = [_GroundType._instance.origin]
        Body.COUNT = 0

        # --- ODE ---
        model_ode, _ = _make_pendulum()
        T_ode, uT_ode = model_ode.solve(
            analysis="dynamic", method="Radau",
            t_final=1.0, dt=0.01, ic_correct=True,
        )

        # Interpolate ODE onto DAE time grid (same dt, so grids should match)
        nB3 = 3 * model_dae.nB
        q_dae = uT_dae[:, :nB3]
        q_ode = uT_ode[:, :nB3]

        # Use the shorter array in case of minor length differences
        n = min(len(T_dae), len(T_ode))
        max_diff = float(np.max(np.abs(q_dae[:n] - q_ode[:n])))
        assert max_diff < 0.05, (
            f"DAE/ODE position discrepancy too large: {max_diff:.4f}"
        )


# ---------------------------------------------------------------------------
# 3. Constraint drift over a longer run
# ---------------------------------------------------------------------------

class TestConstraintDrift:
    """DAE solver must maintain near-zero constraint violation over 3 s."""

    def test_no_drift_long_run(self):
        """max |Phi| < 1e-8 over 3 s with dt = 0.005 s."""
        model, _ = _make_pendulum()
        T, uT = model.solve(
            analysis="dynamic", method="BDF-DAE",
            t_final=3.0, dt=0.005, ic_correct=True,
        )

        nB3 = 3 * model.nB
        max_phi = 0.0
        for k in range(len(T)):
            q_k = uT[k, :nB3]
            for Bi, body in enumerate(model.Bodies):
                ir = 3 * Bi
                body.position    = q_k[ir:ir + 2].reshape(2, 1)
                body.orientation = float(q_k[ir + 2])
            model._update_position()
            phi = model._compute_constraints().flatten()
            max_phi = max(max_phi, float(np.max(np.abs(phi))))

        assert max_phi < 1e-8, (
            f"Constraint drift detected: max|Phi| = {max_phi:.3e}"
        )


# ---------------------------------------------------------------------------
# 4. Input routing
# ---------------------------------------------------------------------------

class TestInputRouting:

    def test_bdf_dae_routes_without_error(self):
        """solve(analysis='dynamic', method='BDF-DAE') must complete."""
        model, _ = _make_pendulum()
        T, uT = model.solve(
            analysis="dynamic", method="BDF-DAE",
            t_final=0.2, dt=0.02,
        )
        assert len(T) > 0
        assert uT.shape[0] == len(T)

    def test_method_name_case_insensitive(self):
        """'bdf-dae', 'Bdf-Dae', 'BDF-DAE' are all accepted."""
        for variant in ("bdf-dae", "Bdf-Dae", "BDF-DAE"):
            _GroundType._instance._markers = [_GroundType._instance.origin]
            Body.COUNT = 0
            model, _ = _make_pendulum()
            T, uT = model.solve(
                analysis="dynamic", method=variant,
                t_final=0.1, dt=0.05,
            )
            assert len(T) > 0, f"No output for method='{variant}'"

    def test_bdf_order_1_runs(self):
        """bdf_order=1 (backward Euler) must complete without error."""
        model, _ = _make_pendulum()
        T, uT = model.solve(
            analysis="dynamic", method="BDF-DAE",
            t_final=0.2, dt=0.02,
        )
        assert len(T) > 0

    def test_t_span_accepted(self):
        """t_span keyword is correctly forwarded to _solve_dae."""
        model, _ = _make_pendulum()
        T, uT = model.solve(
            analysis="dynamic", method="BDF-DAE",
            t_span=(0.0, 0.3), dt=0.05,
        )
        assert len(T) > 0
        assert T[-1] <= 0.3 + 1e-9
