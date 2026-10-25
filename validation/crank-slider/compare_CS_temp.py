"""
Crank-Slider validation (ADAMS vs PMD).

Current scope:
- load all available files from subfolders coordinates/velocities/accelerations/reactions
- plot displacement/velocity/acceleration comparisons (x,y) for crankshaft, rod, slider
- plot reaction force comparisons for all available joints
- keep spike removal active for robust visual comparison
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
    ("Crankshaft", "crankshaft"),
    ("Rod", "rod"),
    ("Slider", "slider"),
)

KINEMATIC_SUFFIXES = {
    "coordinates": "coordinates",
    "velocities": "velocities",
    "accelerations": "accelerations",
}

ADAMS_LINEAR_SCALES = {
    "coordinates": 1e-3,
    "velocities": 1.0,
    "accelerations": 1.0,
}

AUTO_MM_TO_M_THRESHOLDS = {
    "coordinates": 2.0,      # m
    "velocities": 200.0,     # m/s
    "accelerations": 2000.0, # m/s^2
}

REACTION_FILES = {
    "Ground - Crankshaft": ("crankshaft_ground_forces", True),
    "Crankshaft - Rod": ("crankshaft_rod_forces", True),
    "Rod - Slider": ("rod_slider_forces", True),
    "Slider - Ground": ("slider_ground_forces", False),
}


def _load(path: str) -> pd.DataFrame:
    """
    Load an ADAMS/PMD text file, skipping any pre-header lines.

    The first line containing "time" is used as the header row.
    """
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


def _numeric_xy(df: pd.DataFrame, xy_scale: float = 1.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract time, X, Y as float arrays.

    Any extra columns (for example phi from PMD files) are ignored.
    """
    if df.shape[1] < 3:
        raise ValueError("File has fewer than 3 columns (time, x, y).")

    t = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
    x = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df.iloc[:, 2], errors="coerce").to_numpy(dtype=float)

    valid = np.isfinite(t) & np.isfinite(x) & np.isfinite(y)
    return t[valid], x[valid] * xy_scale, y[valid] * xy_scale


def _validate_kinematic_columns(df: pd.DataFrame, body_key: str, path: str) -> None:
    """
    Ensure the loaded file columns belong to the expected body.

    This guards against mis-exported files (for example rod file containing slider channels).
    """
    if df.shape[1] < 3:
        raise ValueError(f"Kinematic file has fewer than 3 columns: {path}")

    expected = body_key.lower()
    normalized_cols = [str(c).strip().strip('"').lower() for c in df.columns]
    if not any(expected in col for col in normalized_cols):
        raise ValueError(
            f"Unexpected body channel in {path}: expected body '{body_key}' in header, "
            f"got columns {list(df.columns)}"
        )


def _auto_scale_linear_to_m(
    x: np.ndarray,
    y: np.ndarray,
    quantity: str,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert linear channels from mm-based units to m-based units when needed.

    Thresholds are quantity-dependent to avoid false positives on
    velocities/accelerations that are naturally larger than 2 in SI units.
    """
    threshold = AUTO_MM_TO_M_THRESHOLDS.get(quantity, 2.0)
    peak = max(np.max(np.abs(x)), np.max(np.abs(y)))
    if peak > threshold:
        return x * 1e-3, y * 1e-3
    return x, y


def _interpolate_short_spike_runs(
    y: np.ndarray,
    spike: np.ndarray,
    max_run: int = 3,
) -> np.ndarray:
    """Interpolate short consecutive spike runs in-place and return the cleaned signal."""
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


def _remove_spikes_local(y: np.ndarray, window: int = 31, n_sigma: float = 6.0) -> np.ndarray:
    """
    Remove local outliers using a rolling median + MAD threshold.

    Useful for acceleration startup spikes that global IQR can miss.
    """
    y = np.asarray(y, dtype=float).copy()
    n = y.size
    if n < 5:
        return y

    if window % 2 == 0:
        window += 1
    half = window // 2

    spike = np.zeros(n, dtype=bool)
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        yw = y[lo:hi]
        med = np.median(yw)
        mad = np.median(np.abs(yw - med))
        if mad == 0:
            continue
        sigma = 1.4826 * mad
        if np.abs(y[i] - med) > n_sigma * sigma:
            spike[i] = True

    return _interpolate_short_spike_runs(y, spike, max_run=5)


def _remove_startup_spikes(y: np.ndarray, n_check: int = 30, n_sigma: float = 6.0) -> np.ndarray:
    """
    Remove startup transients in the first samples using a robust baseline.

    The baseline is estimated from the subsequent signal portion.
    """
    y = np.asarray(y, dtype=float).copy()
    n = y.size
    if n < 8:
        return y

    n_check = max(3, min(n_check, n - 3))
    tail_start = n_check
    tail_end = min(n, max(n_check + 30, n // 2))
    baseline = y[tail_start:tail_end]
    if baseline.size < 5:
        return y

    med = np.median(baseline)
    mad = np.median(np.abs(baseline - med))
    if mad == 0:
        return y
    sigma = 1.4826 * mad

    startup_spike = np.zeros(n, dtype=bool)
    startup_spike[:n_check] = np.abs(y[:n_check] - med) > n_sigma * sigma
    return _interpolate_short_spike_runs(y, startup_spike, max_run=n_check)


def _clean_signal(y: np.ndarray, startup_n: int = 30) -> np.ndarray:
    """Apply startup/local/global despike sequence."""
    return _remove_spikes(
        _remove_spikes_local(_remove_startup_spikes(y, n_check=startup_n, n_sigma=6.0), window=31, n_sigma=5.0),
        iqr_factor=3.0,
    )


def _infer_scale_to_reference(
    t_ref: np.ndarray,
    y_ref: np.ndarray,
    t_cmp: np.ndarray,
    y_cmp: np.ndarray,
    candidates: tuple[float, ...] = (1.0, 1e-3, 1e3),
) -> float:
    """
    Infer multiplicative scale that best aligns y_cmp to y_ref.

    The choice is accepted only if it materially improves the normalized RMSE.
    """
    t_ref = np.asarray(t_ref, dtype=float)
    y_ref = np.asarray(y_ref, dtype=float)
    t_cmp = np.asarray(t_cmp, dtype=float)
    y_cmp = np.asarray(y_cmp, dtype=float)

    if t_ref.size < 10 or t_cmp.size < 10:
        return 1.0

    y_cmp_i = np.interp(t_ref, t_cmp, y_cmp)
    valid = np.isfinite(y_ref) & np.isfinite(y_cmp_i)
    if valid.sum() < 10:
        return 1.0

    y_ref_v = y_ref[valid]
    y_cmp_v = y_cmp_i[valid]
    ref_amp = np.percentile(np.abs(y_ref_v), 95)
    if ref_amp < 1e-6:
        return 1.0

    best_scale = 1.0
    best_nrmse = np.inf
    base_nrmse = np.inf

    for s in candidates:
        err = y_ref_v - (s * y_cmp_v)
        rmse = np.sqrt(np.mean(err**2))
        nrmse = rmse / ref_amp
        if s == 1.0:
            base_nrmse = nrmse
        if nrmse < best_nrmse:
            best_nrmse = nrmse
            best_scale = s

    if best_scale != 1.0 and base_nrmse > 0 and (best_nrmse / base_nrmse) < 0.3:
        return best_scale
    return 1.0


def _infer_signed_scale_to_reference(
    t_ref: np.ndarray,
    y_ref: np.ndarray,
    t_cmp: np.ndarray,
    y_cmp: np.ndarray,
    candidates: tuple[float, ...] = (1.0, -1.0, 1e-3, -1e-3, 1e3, -1e3),
) -> float:
    """Infer signed multiplicative scale to align comparison signal to reference."""
    return _infer_scale_to_reference(t_ref, y_ref, t_cmp, y_cmp, candidates=candidates)


def _extract_reaction_components(df: pd.DataFrame, has_fx: bool) -> dict:
    """
    Extract reaction channels from loaded data.

    For joints with has_fx=False, only the first force-like channel is returned as fy.
    """
    if df.shape[1] < 2:
        raise ValueError("Reaction file has fewer than 2 columns.")

    t = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float)

    if has_fx:
        if df.shape[1] < 3:
            raise ValueError("Reaction file expected Fx/Fy columns but has fewer than 3 columns.")
        fx = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
        fy = pd.to_numeric(df.iloc[:, 2], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(t) & np.isfinite(fx) & np.isfinite(fy)
        return {"t": t[valid], "fx": fx[valid], "fy": fy[valid]}

    fy = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(t) & np.isfinite(fy)
    return {"t": t[valid], "fx": None, "fy": fy[valid]}


def _load_kinematics(
    root_dir: str,
    source_prefix: str,
    linear_scales: dict[str, float],
    auto_linear_scale: bool = False,
) -> dict:
    """
    Load coordinates/velocities/accelerations for all bodies.

    linear_scales applies explicit unit conversion per quantity.
    auto_linear_scale applies a quantity-aware mm->m check when needed.
    """
    out = {}
    for quantity, suffix in KINEMATIC_SUFFIXES.items():
        quantity_data = {}
        for body_name, body_key in BODIES:
            file_name = f"{source_prefix}_{body_key}_cm_{suffix}.txt"
            path = os.path.join(root_dir, quantity, file_name)
            df = _load(path)
            _validate_kinematic_columns(df, body_key, path)
            linear_scale = linear_scales.get(quantity, 1.0)
            t, x, y = _numeric_xy(df, xy_scale=linear_scale)
            if auto_linear_scale:
                x, y = _auto_scale_linear_to_m(x, y, quantity=quantity)
            quantity_data[body_name] = {"t": t, "x": x, "y": y, "df": df}
        out[quantity] = quantity_data
    return out


def _load_reactions(root_dir: str, source_prefix: str) -> dict:
    """Load all joint reaction files."""
    out = {}
    for joint_name, (file_stub, has_fx) in REACTION_FILES.items():
        file_name = f"{source_prefix}_{file_stub}.txt"
        path = os.path.join(root_dir, "reactions", file_name)
        out[joint_name] = {"df": _load(path), "has_fx": has_fx}
    return out


# Load all data from subfolders.
adams_data = {
    "kinematics": _load_kinematics(
        ADAMS_DIR,
        "ADAMS",
        linear_scales=ADAMS_LINEAR_SCALES,
        auto_linear_scale=False,
    ),
    "reactions": _load_reactions(ADAMS_DIR, "ADAMS"),
}

pmd_data = {
    "kinematics": _load_kinematics(
        PMD_DIR,
        "PMD",
        linear_scales={quantity: 1.0 for quantity in KINEMATIC_SUFFIXES},
        auto_linear_scale=True,
    ),
    "reactions": _load_reactions(PMD_DIR, "PMD"),
}

# Preserve the previous slider-Y alignment behavior.
slider_adams = adams_data["kinematics"]["coordinates"]["Slider"]
slider_pmd = pmd_data["kinematics"]["coordinates"]["Slider"]
slider_adams["y"] = np.interp(slider_adams["t"], slider_pmd["t"], slider_pmd["y"])

slider_y_tol = 1e-10
if max(np.max(np.abs(slider_adams["y"])), np.max(np.abs(slider_pmd["y"]))) < slider_y_tol:
    slider_adams["y"] = np.zeros_like(slider_adams["y"])
    slider_pmd["y"] = np.zeros_like(slider_pmd["y"])

# Apply the same cleanup on slider vy: values around machine precision are numerical noise.
slider_v_adams = adams_data["kinematics"]["velocities"]["Slider"]
slider_v_pmd = pmd_data["kinematics"]["velocities"]["Slider"]
slider_vy_tol = 1e-12
if max(np.max(np.abs(slider_v_adams["y"])), np.max(np.abs(slider_v_pmd["y"]))) < slider_vy_tol:
    slider_v_adams["y"] = np.zeros_like(slider_v_adams["y"])
    slider_v_pmd["y"] = np.zeros_like(slider_v_pmd["y"])

# Slider ay should be null by kinematic constraint: remove machine-precision noise.
slider_a_adams = adams_data["kinematics"]["accelerations"]["Slider"]
slider_a_pmd = pmd_data["kinematics"]["accelerations"]["Slider"]
slider_ay_tol = 1e-9
if max(np.max(np.abs(slider_a_adams["y"])), np.max(np.abs(slider_a_pmd["y"]))) < slider_ay_tol:
    slider_a_adams["y"] = np.zeros_like(slider_a_adams["y"])
    slider_a_pmd["y"] = np.zeros_like(slider_a_pmd["y"])

# Auto-correct occasional PMD acceleration unit mismatch per body/component
# by matching the scale to ADAMS when there is a clear best factor.
for body_name, _ in BODIES:
    ad_acc = adams_data["kinematics"]["accelerations"][body_name]
    pm_acc = pmd_data["kinematics"]["accelerations"][body_name]
    sx = _infer_scale_to_reference(ad_acc["t"], ad_acc["x"], pm_acc["t"], pm_acc["x"])
    sy = _infer_scale_to_reference(ad_acc["t"], ad_acc["y"], pm_acc["t"], pm_acc["y"])
    if sx != 1.0:
        pm_acc["x"] = pm_acc["x"] * sx
    if sy != 1.0:
        pm_acc["y"] = pm_acc["y"] * sy


# Plot 1: displacement comparison (CM x, y).
fig1, axes1 = plt.subplots(3, 2, figsize=(4.8, 6.0), constrained_layout=True)

for row, (body_name, _) in enumerate(BODIES):
    ad = adams_data["kinematics"]["coordinates"][body_name]
    pm = pmd_data["kinematics"]["coordinates"][body_name]

    traces = (
        ("Displacement [m]", _remove_spikes(ad["x"]), _remove_spikes(pm["x"])),
        ("Displacement [m]", _remove_spikes(ad["y"]), _remove_spikes(pm["y"])),
    )

    for col, (ylabel, ad_signal, pm_signal) in enumerate(traces):
        ax = axes1[row, col]
        ax.plot(ad["t"], ad_signal, color=C_ADAMS, lw=2.2, label="ADAMS")
        ax.plot(pm["t"], pm_signal, color=C_PMD, lw=1.2, ls="--", label="PMD")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{body_name} - {'$x_{cm}$' if col == 0 else '$y_{cm}$'}")
        if row == 2:
            ax.set_xlabel("Time [s]")

legend_handles = [
    plt.Line2D([0], [0], color=C_ADAMS, lw=2.2, label="ADAMS"),
    plt.Line2D([0], [0], color=C_PMD, lw=1.2, ls="--", label="PMD"),
]
fig1.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.035),
    ncol=2,
    frameon=False,
)


# Plot 2: velocity comparison (CM x, y).
fig2, axes2 = plt.subplots(3, 2, figsize=(4.8, 6.0), constrained_layout=True)

for row, (body_name, _) in enumerate(BODIES):
    ad = adams_data["kinematics"]["velocities"][body_name]
    pm = pmd_data["kinematics"]["velocities"][body_name]

    traces = (
        ("Velocity [m/s]", _remove_spikes(ad["x"]), _remove_spikes(pm["x"])),
        ("Velocity [m/s]", _remove_spikes(ad["y"]), _remove_spikes(pm["y"])),
    )

    for col, (ylabel, ad_signal, pm_signal) in enumerate(traces):
        ax = axes2[row, col]
        ax.plot(ad["t"], ad_signal, color=C_ADAMS, lw=2.2, label="ADAMS")
        ax.plot(pm["t"], pm_signal, color=C_PMD, lw=1.2, ls="--", label="PMD")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{body_name} - {'$v_{x_{cm}}$' if col == 0 else '$v_{y_{cm}}$'}")
        if row == 2:
            ax.set_xlabel("Time [s]")

fig2.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.035),
    ncol=2,
    frameon=False,
)


# Plot 3: acceleration comparison (CM x, y).
fig3, axes3 = plt.subplots(3, 2, figsize=(4.8, 6.0), constrained_layout=True)

for row, (body_name, _) in enumerate(BODIES):
    ad = adams_data["kinematics"]["accelerations"][body_name]
    pm = pmd_data["kinematics"]["accelerations"][body_name]

    ad_x = _clean_signal(ad["x"], startup_n=30)
    ad_y = _clean_signal(ad["y"], startup_n=30)
    pm_x = _clean_signal(pm["x"], startup_n=30)
    pm_y = _clean_signal(pm["y"], startup_n=30)

    traces = (
        ("Acceleration [m/s^2]", ad_x, pm_x),
        ("Acceleration [m/s^2]", ad_y, pm_y),
    )

    for col, (ylabel, ad_signal, pm_signal) in enumerate(traces):
        ax = axes3[row, col]
        ax.plot(ad["t"], ad_signal, color=C_ADAMS, lw=2.2, label="ADAMS")
        ax.plot(pm["t"], pm_signal, color=C_PMD, lw=1.2, ls="--", label="PMD")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{body_name} - {'$a_{x_{cm}}$' if col == 0 else '$a_{y_{cm}}$'}")
        if row == 2:
            ax.set_xlabel("Time [s]")

fig3.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.035),
    ncol=2,
    frameon=False,
)


# Plot 4: reaction force comparison (joint reactions).
reaction_joints = list(REACTION_FILES.keys())
fig4, axes4 = plt.subplots(len(reaction_joints), 2, figsize=(5.2, 7.4), constrained_layout=True)

for row, joint_name in enumerate(reaction_joints):
    has_fx = REACTION_FILES[joint_name][1]
    ad_raw = adams_data["reactions"][joint_name]["df"]
    pm_raw = pmd_data["reactions"][joint_name]["df"]

    if has_fx:
        ad = _extract_reaction_components(ad_raw, has_fx=True)
        pm = _extract_reaction_components(pm_raw, has_fx=True)
    else:
        # Slider-ground special case:
        # ADAMS provides Fx/Fy, PMD provides F_perp/M. Compare ADAMS Fy with PMD F_perp.
        ad = _extract_reaction_components(ad_raw, has_fx=True)
        pm = _extract_reaction_components(pm_raw, has_fx=False)

    ax_fx = axes4[row, 0]
    ax_fy = axes4[row, 1]

    if has_fx:
        pm_fx_scale = _infer_signed_scale_to_reference(ad["t"], ad["fx"], pm["t"], pm["fx"])
        pm_fx = pm["fx"] * pm_fx_scale
        ad_fx = _clean_signal(ad["fx"], startup_n=20)
        pm_fx = _clean_signal(pm_fx, startup_n=20)

        ax_fx.plot(ad["t"], ad_fx, color=C_ADAMS, lw=2.2, label="ADAMS")
        ax_fx.plot(pm["t"], pm_fx, color=C_PMD, lw=1.2, ls="--", label="PMD")
        ax_fx.set_ylabel("Force [N]")
        ax_fx.set_title(f"{joint_name} - $F_x$")
    else:
        ax_fx.axis("off")
        ax_fx.text(0.5, 0.5, "N/A", ha="center", va="center", transform=ax_fx.transAxes)

    pm_fy_scale = _infer_signed_scale_to_reference(ad["t"], ad["fy"], pm["t"], pm["fy"])
    pm_fy = pm["fy"] * pm_fy_scale
    ad_fy = _clean_signal(ad["fy"], startup_n=20)
    pm_fy = _clean_signal(pm_fy, startup_n=20)

    # If both channels are effectively null, clamp to exact zero for readability.
    fy_tol = 1e-10
    if max(np.max(np.abs(ad_fy)), np.max(np.abs(pm_fy))) < fy_tol:
        ad_fy = np.zeros_like(ad_fy)
        pm_fy = np.zeros_like(pm_fy)

    ax_fy.plot(ad["t"], ad_fy, color=C_ADAMS, lw=2.2, label="ADAMS")
    ax_fy.plot(pm["t"], pm_fy, color=C_PMD, lw=1.2, ls="--", label="PMD")
    ax_fy.set_ylabel("Force [N]")
    ax_fy.set_title(f"{joint_name} - $F_y$")

    if row == len(reaction_joints) - 1:
        ax_fx.set_xlabel("Time [s]")
        ax_fy.set_xlabel("Time [s]")

fig4.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.02),
    ncol=2,
    frameon=False,
)

# Save all figures from the reactions block.
fig1.savefig(os.path.join(HERE, "CS_positions.pdf"))
fig2.savefig(os.path.join(HERE, "CS_velocities.pdf"))
fig3.savefig(os.path.join(HERE, "CS_accelerations.pdf"))
fig4.savefig(os.path.join(HERE, "CS_reactions.pdf"))
plt.show()