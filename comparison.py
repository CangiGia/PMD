"""
Compare MATLAB vs Python simulation results for PMD models.

Reads result files from:
  - MATLAB: _references/Nikravesh - .../DAP_BC_12_2018/_results/*.txt
  - Python: examples/results/_test_*.txt

Produces:
  - Console summary with max absolute/relative errors per model
  - Optional plots comparing time histories
  
Usage:
    python compare_results.py
    python compare_results.py --plot
    python compare_results.py --models AA SP CB
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR = Path(__file__).parent

MATLAB_DIR = (
    BASE_DIR
    / "_references"
    / "Nikravesh - Planar Multibody Dynamics 2nd Edition"
    / "DAP_BC_12_2018"
    / "_results"
)

PYTHON_DIR = BASE_DIR / "examples" / "results"

# =============================================================================
# MODEL NAME MAPPING
# =============================================================================

# MATLAB filename -> Python filename
MODEL_MAP = {
    "AA":     "_test_AA",
    "Cart_A": "_test_Cart_A",
    "Cart_B": "_test_Cart_B",
    "Cart_C": "_test_Cart_C",
    "Cart_D": "_test_Cart_D",
    "CB":     "_test_CB",
    "MP_A":   "_test_MP_A",
    "MP_B":   "_test_MP_B",
    "MP_C":   "_test_MP_C",
    "Rod":    "_test_Rod",
    "SP":     "_test_SP",
}


# =============================================================================
# FILE READING
# =============================================================================

def read_result_file(filepath: Path) -> tuple[np.ndarray, list[str]]:
    """
    Read a result file. Auto-detects delimiter (comma, tab, or space).
    
    Returns:
        data: numpy array (n_rows, n_cols)
        headers: list of column names (empty if no header)
    """
    headers = []
    skip_rows = 0

    with open(filepath, "r") as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()

    # Detect delimiter from the first data line
    test_line = first_line
    try:
        [float(x) for x in first_line.split()]
    except ValueError:
        # Check if it's comma-separated numbers
        try:
            [float(x) for x in first_line.split(",")]
            test_line = first_line  # first line IS data, comma-separated
        except ValueError:
            # It's a header line
            headers = first_line.replace(",", "\t").split()
            skip_rows = 1
            test_line = second_line

    # Detect delimiter from data line
    if "," in test_line:
        delimiter = ","
    elif "\t" in test_line:
        delimiter = "\t"
    else:
        delimiter = None  # whitespace (numpy default)

    data = np.loadtxt(filepath, skiprows=skip_rows, delimiter=delimiter)
    return data, headers

# =============================================================================
# COMPARISON LOGIC
# =============================================================================

def compare_model(
    model_name: str,
    matlab_data: np.ndarray,
    python_data: np.ndarray,
    matlab_headers: list[str],
    python_headers: list[str],
) -> dict:
    """
    Compare MATLAB and Python results for a single model.
    
    Handles:
      - Different number of time steps (interpolates to common time base)
      - Different number of columns (compares only common columns)
      - Different time ranges (compares only overlapping range)
    
    Returns:
        Dictionary with comparison metrics.
    """
    results = {"model": model_name, "status": "OK", "details": []}

    # --- Extract time vectors ---
    t_matlab = matlab_data[:, 0]
    t_python = python_data[:, 0]

    # --- Find overlapping time range ---
    t_start = max(t_matlab[0], t_python[0])
    t_end = min(t_matlab[-1], t_python[-1])

    if t_start >= t_end:
        results["status"] = "ERROR"
        results["details"].append("No overlapping time range!")
        return results

    results["t_range_matlab"] = (t_matlab[0], t_matlab[-1])
    results["t_range_python"] = (t_python[0], t_python[-1])
    results["t_range_common"] = (t_start, t_end)

    # --- Common time base (use finer grid) ---
    dt_matlab = np.median(np.diff(t_matlab))
    dt_python = np.median(np.diff(t_python))
    dt = min(dt_matlab, dt_python)
    t_common = np.arange(t_start, t_end + dt / 2, dt)

    # --- Number of data columns (excluding time) ---
    n_cols_matlab = matlab_data.shape[1] - 1
    n_cols_python = python_data.shape[1] - 1
    n_cols_compare = min(n_cols_matlab, n_cols_python)

    if n_cols_matlab != n_cols_python:
        results["details"].append(
            f"Column count mismatch: MATLAB={n_cols_matlab+1}, Python={n_cols_python+1}. "
            f"Comparing first {n_cols_compare} data columns."
        )

    # --- Interpolate both to common time base ---
    from scipy.interpolate import interp1d

    matlab_interp = np.zeros((len(t_common), n_cols_compare))
    python_interp = np.zeros((len(t_common), n_cols_compare))

    for col in range(n_cols_compare):
        f_m = interp1d(t_matlab, matlab_data[:, col + 1], kind="linear",
                       fill_value="extrapolate")
        f_p = interp1d(t_python, python_data[:, col + 1], kind="linear",
                       fill_value="extrapolate")
        matlab_interp[:, col] = f_m(t_common)
        python_interp[:, col] = f_p(t_common)

    # --- Build column name list (needed before error computation) ---
    col_names = []
    for col in range(n_cols_compare):
        if python_headers and col + 1 < len(python_headers):
            col_names.append(python_headers[col + 1])
        elif matlab_headers and col + 1 < len(matlab_headers):
            col_names.append(matlab_headers[col + 1])
        else:
            col_names.append(f"col_{col+1}")

    # --- Compute errors (angle-aware for phi columns) ---
    # For columns whose name ends with '_p', normalize the difference to [-π, π]
    # to remove 2π-wrapping artifacts (e.g. MATLAB stores -0.037 rad while Python
    # stores the equivalent +6.246 rad after initial-condition unwrapping).
    signed_diff = matlab_interp - python_interp
    for col in range(n_cols_compare):
        if col_names[col].endswith("_p"):
            signed_diff[:, col] = (
                (signed_diff[:, col] + np.pi) % (2 * np.pi) - np.pi
            )
    abs_error = np.abs(signed_diff)

    # Relative error (avoid division by zero)
    scale = np.maximum(np.abs(matlab_interp), 1e-10)
    rel_error = abs_error / scale

    # --- Per-column metrics ---
    col_metrics = []
    for col in range(n_cols_compare):
        col_name = col_names[col]

        max_abs = np.max(abs_error[:, col])
        mean_abs = np.mean(abs_error[:, col])
        max_rel = np.max(rel_error[:, col])
        mean_rel = np.mean(rel_error[:, col])
        rms = np.sqrt(np.mean(signed_diff[:, col] ** 2))

        col_metrics.append({
            "name": col_name,
            "max_abs_error": max_abs,
            "mean_abs_error": mean_abs,
            "max_rel_error": max_rel,
            "mean_rel_error": mean_rel,
            "rms_error": rms,
        })

    results["n_cols_compared"] = n_cols_compare
    results["n_timesteps"] = len(t_common)
    results["dt"] = dt
    results["columns"] = col_metrics
    results["max_abs_error_global"] = np.max(abs_error)
    results["max_rel_error_global"] = np.max(rel_error)
    results["rms_error_global"] = np.sqrt(np.mean(signed_diff ** 2))

    # Store for plotting
    results["t_common"] = t_common
    results["matlab_interp"] = matlab_interp
    results["python_interp"] = python_interp

    # --- Status classification ---
    if results["max_abs_error_global"] < 1e-6:
        results["status"] = "EXCELLENT"
    elif results["max_abs_error_global"] < 1e-3:
        results["status"] = "GOOD"
    elif results["max_abs_error_global"] < 1e-1:
        results["status"] = "ACCEPTABLE"
    else:
        results["status"] = "POOR"

    return results


# =============================================================================
# REPORTING
# =============================================================================

STATUS_ICONS = {
    "EXCELLENT": "🟢",
    "GOOD":      "🟡",
    "ACCEPTABLE": "🟠",
    "POOR":      "🔴",
    "ERROR":     "❌",
    "MISSING":   "⚪",
}


def print_summary(all_results: list[dict]):
    """Print a summary table of all model comparisons."""
    
    print("\n" + "=" * 90)
    print("  PMD — MATLAB vs Python Comparison Summary")
    print("=" * 90)

    # --- Summary table ---
    print(f"\n{'Model':<12} {'Status':<14} {'Max |Err|':<14} {'RMS Err':<14} {'Max Rel Err':<14} {'Cols':<6}")
    print("-" * 76)

    for r in all_results:
        icon = STATUS_ICONS.get(r["status"], "?")
        if r["status"] in ("ERROR", "MISSING"):
            print(f"{r['model']:<12} {icon} {r['status']:<11} {'—':<14} {'—':<14} {'—':<14} {'—':<6}")
            for d in r.get("details", []):
                print(f"{'':>14} ↳ {d}")
        else:
            print(
                f"{r['model']:<12} {icon} {r['status']:<11} "
                f"{r['max_abs_error_global']:<14.2e} "
                f"{r['rms_error_global']:<14.2e} "
                f"{r['max_rel_error_global']:<14.2e} "
                f"{r['n_cols_compared']:<6}"
            )

    # --- Detailed per-column breakdown ---
    print("\n" + "=" * 90)
    print("  Detailed Per-Column Errors")
    print("=" * 90)

    for r in all_results:
        if r["status"] in ("ERROR", "MISSING"):
            continue

        print(f"\n┌─ {r['model']} ({r['n_timesteps']} timesteps, dt={r['dt']:.4f}s)")
        if r.get("details"):
            for d in r["details"]:
                print(f"│  ⚠ {d}")

        print(f"│  {'Column':<16} {'Max |Err|':<14} {'Mean |Err|':<14} {'RMS':<14} {'Max Rel':<14}")
        print(f"│  {'-'*70}")

        for cm in r["columns"]:
            print(
                f"│  {cm['name']:<16} "
                f"{cm['max_abs_error']:<14.2e} "
                f"{cm['mean_abs_error']:<14.2e} "
                f"{cm['rms_error']:<14.2e} "
                f"{cm['max_rel_error']:<14.2e}"
            )
        print("└" + "─" * 74)

    # --- Overall verdict ---
    statuses = [r["status"] for r in all_results]
    n_total = len(all_results)
    n_excellent = statuses.count("EXCELLENT")
    n_good = statuses.count("GOOD")
    n_acceptable = statuses.count("ACCEPTABLE")
    n_poor = statuses.count("POOR")
    n_error = statuses.count("ERROR")
    n_missing = statuses.count("MISSING")

    print(f"\n{'='*90}")
    print(f"  OVERALL: {n_total} models | "
          f"🟢 {n_excellent} excellent | 🟡 {n_good} good | "
          f"🟠 {n_acceptable} acceptable | 🔴 {n_poor} poor | "
          f"❌ {n_error} errors | ⚪ {n_missing} missing")
    print(f"{'='*90}\n")


# =============================================================================
# PLOTTING
# =============================================================================

def plot_comparison(result: dict):
    """Plot MATLAB vs Python time histories for a single model."""
    
    if result["status"] in ("ERROR", "MISSING"):
        print(f"  Skipping plot for {result['model']} ({result['status']})")
        return

    t = result["t_common"]
    matlab = result["matlab_interp"]
    python = result["python_interp"]
    n_cols = result["n_cols_compared"]

    # Group columns by body (3 columns per body: x, y, p)
    n_bodies = n_cols // 3 if n_cols >= 3 else 1
    cols_per_body = 3 if n_cols >= 3 else n_cols

    fig, axes = plt.subplots(
        n_bodies, cols_per_body,
        figsize=(5 * cols_per_body, 3 * n_bodies),
        squeeze=False,
        sharex=True,
    )
    fig.suptitle(f"Model: {result['model']} — MATLAB vs Python", fontsize=14, fontweight="bold")

    labels = ["x [m]", "y [m]", "φ [rad]"]

    for body_idx in range(n_bodies):
        for sub_idx in range(cols_per_body):
            col = body_idx * cols_per_body + sub_idx
            if col >= n_cols:
                break

            ax = axes[body_idx, sub_idx]

            # Get column name
            col_name = result["columns"][col]["name"]

            ax.plot(t, matlab[:, col], "b-", linewidth=1.2, label="MATLAB", alpha=0.8)
            ax.plot(t, python[:, col], "r--", linewidth=1.2, label="Python", alpha=0.8)

            ax.set_ylabel(col_name if col_name else labels[sub_idx % 3])
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)

            if body_idx == 0 and sub_idx == 0:
                ax.set_title(f"Body {body_idx + 1}")
            elif sub_idx == 0:
                ax.set_title(f"Body {body_idx + 1}")

    for ax in axes[-1, :]:
        ax.set_xlabel("Time [s]")

    plt.tight_layout()

    # --- Error plot ---
    fig2, axes2 = plt.subplots(
        n_bodies, cols_per_body,
        figsize=(5 * cols_per_body, 2.5 * n_bodies),
        squeeze=False,
        sharex=True,
    )
    fig2.suptitle(
        f"Model: {result['model']} — Absolute Error |MATLAB - Python|",
        fontsize=14, fontweight="bold",
    )

    for body_idx in range(n_bodies):
        for sub_idx in range(cols_per_body):
            col = body_idx * cols_per_body + sub_idx
            if col >= n_cols:
                break

            ax = axes2[body_idx, sub_idx]
            col_name = result["columns"][col]["name"]
            error = np.abs(matlab[:, col] - python[:, col])

            ax.semilogy(t, error + 1e-16, "k-", linewidth=0.8)
            ax.set_ylabel(f"|Δ{col_name}|" if col_name else f"|Δ{labels[sub_idx % 3]}|")
            ax.grid(True, alpha=0.3)

            if sub_idx == 0:
                ax.set_title(f"Body {body_idx + 1}")

    for ax in axes2[-1, :]:
        ax.set_xlabel("Time [s]")

    plt.tight_layout()


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Compare MATLAB vs Python PMD simulation results."
    )
    parser.add_argument(
        "--plot", action="store_true",
        help="Show comparison plots for each model."
    )
    parser.add_argument(
        "--models", nargs="*", default=None,
        help="List of model names to compare (default: all). E.g.: --models AA SP CB"
    )
    args = parser.parse_args()

    # --- Validate directories ---
    if not MATLAB_DIR.exists():
        print(f"ERROR: MATLAB results directory not found:\n  {MATLAB_DIR}")
        sys.exit(1)
    if not PYTHON_DIR.exists():
        print(f"ERROR: Python results directory not found:\n  {PYTHON_DIR}")
        sys.exit(1)

    # --- Select models ---
    if args.models:
        models_to_compare = {k: v for k, v in MODEL_MAP.items() if k in args.models}
        unknown = set(args.models) - set(MODEL_MAP.keys())
        if unknown:
            print(f"WARNING: Unknown models ignored: {unknown}")
            print(f"  Available: {list(MODEL_MAP.keys())}")
    else:
        models_to_compare = MODEL_MAP

    # --- Run comparisons ---
    all_results = []

    for matlab_name, python_name in sorted(models_to_compare.items()):
        matlab_file = MATLAB_DIR / f"{matlab_name}.txt"
        python_file = PYTHON_DIR / f"{python_name}.txt"

        print(f"Comparing {matlab_name}...", end=" ")

        # Check file existence
        if not matlab_file.exists():
            result = {"model": matlab_name, "status": "MISSING",
                      "details": [f"MATLAB file not found: {matlab_file}"]}
            all_results.append(result)
            print("MISSING (MATLAB)")
            continue

        if not python_file.exists():
            result = {"model": matlab_name, "status": "MISSING",
                      "details": [f"Python file not found: {python_file}"]}
            all_results.append(result)
            print("MISSING (Python)")
            continue

        try:
            matlab_data, matlab_headers = read_result_file(matlab_file)
            python_data, python_headers = read_result_file(python_file)

            result = compare_model(
                matlab_name, matlab_data, python_data,
                matlab_headers, python_headers,
            )
            all_results.append(result)
            print(f"{STATUS_ICONS.get(result['status'], '?')} {result['status']}")

        except Exception as e:
            result = {"model": matlab_name, "status": "ERROR",
                      "details": [f"Exception: {e}"]}
            all_results.append(result)
            print(f"ERROR: {e}")

    # --- Print summary ---
    print_summary(all_results)

    # --- Plots ---
    if args.plot:
        for r in all_results:
            plot_comparison(r)
        plt.show()


if __name__ == "__main__":
    main()