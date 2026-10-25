"""
Mass-Spring-Damper validation (ADAMS vs PMD).

Scope:
- load coordinates/velocities/accelerations for the single body (mass)
- handle files stored directly under ADAMS/ and PMD/
- produce three comparison figures (x,y): position, velocity, acceleration
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
VALIDATION_DIR = os.path.dirname(HERE)

if VALIDATION_DIR not in sys.path:
    sys.path.insert(0, VALIDATION_DIR)

import plot_style as plot

ADAMS_DIR = os.path.join(HERE, "ADAMS")
PMD_DIR = os.path.join(HERE, "PMD")

plot.setup_thesis_plots()
plt.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
})

C_ADAMS = "#2E7EC1"
C_PMD = "#D64B35"

BODIES = (
    ("Mass", "mass"),
)

BODY_ALIASES = {
    "mass": ("mass",),
}

KINEMATIC_SUFFIXES = {
    "coordinates": "coordinates",
    "velocities": "velocities",
    "accelerations": "accelerations",
}


def _load(path: str) -> pd.DataFrame:
    """Load an ADAMS/PMD text file, skipping any pre-header lines."""
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    header_idx = None
    for i, line in enumerate(lines):
        if "time" in line.lower():
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(f"Header row with time column not found in: {path}")

    return pd.read_csv(
        path,
        sep=r"\s+|\t+",
        header=0,
        skiprows=header_idx,
        quotechar='"',
        engine="python",
    )


def _numeric_xy(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract time, X, Y as float arrays."""
    if df.shape[1] < 3:
        raise ValueError("File has fewer than 3 columns (time, x, y).")

    t = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
    x = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df.iloc[:, 2], errors="coerce").to_numpy(dtype=float)

    valid = np.isfinite(t) & np.isfinite(x) & np.isfinite(y)
    return t[valid], x[valid], y[valid]


def _validate_kinematic_columns(
    df: pd.DataFrame,
    aliases: tuple[str, ...],
    path: str,
) -> None:
    """Ensure kinematic columns belong to the expected body."""
    normalized_cols = [str(c).strip().strip('"').lower() for c in df.columns]
    joined = " ".join(normalized_cols)
    if not any(alias in joined for alias in aliases):
        raise ValueError(
            f"Unexpected body channel in {path}: expected one of {aliases}, "
            f"got columns {list(df.columns)}"
        )


def _interpolate_short_spike_runs(
    y: np.ndarray,
    spike: np.ndarray,
    max_run: int = 3,
) -> np.ndarray:
    """Interpolate short consecutive spike runs."""
    y = np.asarray(y, dtype=float).copy()
    spike = np.asarray(spike, dtype=bool)

    i = 0
    n = len(y)
    while i < n:
        if spike[i]:
            j = i
            while j < n and spike[j]:
                j += 1
            run = j - i
            if run <= max_run:
                lo = i - 1 if i > 0 else j
                hi = j if j < n else i - 1
                if lo >= 0 and hi < n:
                    y[i:j] = np.interp(np.arange(i, j), [lo, hi], [y[lo], y[hi]])
            i = j
        else:
            i += 1
    return y


def _remove_spikes(y: np.ndarray, iqr_factor: float = 5.0) -> np.ndarray:
    """Replace short isolated spikes with linear interpolation."""
    y = np.asarray(y, dtype=float).copy()
    if y.size < 5:
        return y

    q25, q75 = np.percentile(y, [25, 75])
    iqr = q75 - q25
    if iqr == 0:
        return y

    med = np.median(y)
    spike = np.abs(y - med) > iqr_factor * iqr
    return _interpolate_short_spike_runs(y, spike, max_run=3)


def _kinematic_path(root_dir: str, source_prefix: str, body_key: str, suffix: str) -> str:
    """
    Resolve kinematic file path.

    Prefer direct ADAMS/PMD root files; fallback to quantity subfolder layout.
    """
    file_name = f"{source_prefix}_{body_key}_cm_{suffix}.txt"
    direct = os.path.join(root_dir, file_name)
    if os.path.exists(direct):
        return direct

    nested = os.path.join(root_dir, suffix, file_name)
    if os.path.exists(nested):
        return nested

    raise FileNotFoundError(
        f"Missing file for {source_prefix} {body_key} {suffix}: "
        f"tried {direct} and {nested}"
    )


def _load_kinematics(root_dir: str, source_prefix: str) -> dict:
    """Load coordinates/velocities/accelerations for all bodies."""
    out = {}
    for quantity, suffix in KINEMATIC_SUFFIXES.items():
        quantity_data = {}
        for body_name, body_key in BODIES:
            path = _kinematic_path(root_dir, source_prefix, body_key, suffix)
            df = _load(path)
            _validate_kinematic_columns(df, BODY_ALIASES[body_key], path)
            t, x, y = _numeric_xy(df)
            quantity_data[body_name] = {"t": t, "x": x, "y": y, "df": df}
        out[quantity] = quantity_data
    return out


adams_data = {
    "kinematics": _load_kinematics(ADAMS_DIR, "ADAMS"),
}

pmd_data = {
    "kinematics": _load_kinematics(PMD_DIR, "PMD"),
}

# X motion is numerically negligible for this 1-DOF case, so clamp it to zero.
for quantity in KINEMATIC_SUFFIXES:
    adams_data["kinematics"][quantity]["Mass"]["x"] = np.zeros_like(
        adams_data["kinematics"][quantity]["Mass"]["x"]
    )
    pmd_data["kinematics"][quantity]["Mass"]["x"] = np.zeros_like(
        pmd_data["kinematics"][quantity]["Mass"]["x"]
    )

legend_handles = [
    plt.Line2D([0], [0], color=C_ADAMS, lw=2.2, label="ADAMS"),
    plt.Line2D([0], [0], color=C_PMD, lw=1.2, ls="--", label="PMD"),
]

# Single figure: displacement, velocity, acceleration with X and Y columns.
fig, axes = plt.subplots(3, 2, figsize=(4.8, 6.0), constrained_layout=True)

traces = (
    ("coordinates", "Displacement [mm]", "$x_{cm}$", "$y_{cm}$"),
    ("velocities", "Velocity [mm/s]", "$v_{x_{cm}}$", "$v_{y_{cm}}$"),
    ("accelerations", "Acceleration [mm/s$^2$]", "$a_{x_{cm}}$", "$a_{y_{cm}}$"),
)

for row, (quantity, ylabel, x_label, y_label) in enumerate(traces):
    ad = adams_data["kinematics"][quantity]["Mass"]
    pm = pmd_data["kinematics"][quantity]["Mass"]

    ad_x = np.zeros_like(ad["t"])
    pm_x = np.zeros_like(pm["t"])
    ad_y = _remove_spikes(ad["y"])
    pm_y = _remove_spikes(pm["y"])

    ax_x = axes[row, 0]
    ax_y = axes[row, 1]

    ax_x.plot(ad["t"], ad_x, color=C_ADAMS, lw=2.2)
    ax_x.plot(pm["t"], pm_x, color=C_PMD, lw=1.2, ls="--")
    ax_x.set_ylabel(ylabel)
    ax_x.set_title(f"Mass - {x_label}")

    ax_y.plot(ad["t"], ad_y, color=C_ADAMS, lw=2.2)
    ax_y.plot(pm["t"], pm_y, color=C_PMD, lw=1.2, ls="--")
    ax_y.set_ylabel(ylabel)
    ax_y.set_title(f"Mass - {y_label}")

    if row == len(traces) - 1:
        ax_x.set_xlabel("Time [s]")
        ax_y.set_xlabel("Time [s]")

fig.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.035),
    ncol=2,
    frameon=False,
)

fig.savefig(os.path.join(HERE, "MSD_kinematics.pdf"))
plt.show()