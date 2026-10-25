"""
Four-Bar Linkage validation (ADAMS vs PMD).

Current scope:
- load all available files from subfolders coordinates/velocities/accelerations/reactions
- plot displacement/velocity/acceleration comparisons (x,y) for crank, coupler, follower
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
    ("Crankshaft", "crank"),
    ("Coupler", "coupler"),
    ("Follower", "follower"),
)

BODY_ALIASES = {
    "crank": ("crank", "link_1"),
    "coupler": ("coupler", "link_2"),
    "follower": ("follower", "link_3"),
}

KINEMATIC_SUFFIXES = {
    "coordinates": "coordinates",
    "velocities": "velocities",
    "accelerations": "accelerations",
}

ADAMS_LINEAR_SCALES = {
    "coordinates": 1e-3,
    "velocities": 1e-3,
    "accelerations": 1e-3,
}

REACTION_FILES = {
    "Ground - Crank": ("ground_crank_forces", True),
    "Crank - Coupler": ("crank_coupler_forces", True),
    "Coupler - Follower": ("coupler_follower_forces", True),
    "Follower - Ground": ("follower_ground_forces", True),
}


def _load(path: str) -> pd.DataFrame:
    """Load an ADAMS/PMD text file, skipping pre-header lines."""
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
        sep=r"\t+",
        header=0,
        skiprows=header_idx,
        quotechar='"',
        engine="python",
    )


def _numeric_xy(df: pd.DataFrame, xy_scale: float = 1.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract time, X, Y as float arrays."""
    if df.shape[1] < 3:
        raise ValueError("File has fewer than 3 columns (time, x, y).")

    t = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
    x = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df.iloc[:, 2], errors="coerce").to_numpy(dtype=float)

    valid = np.isfinite(t) & np.isfinite(x) & np.isfinite(y)
    return t[valid], x[valid] * xy_scale, y[valid] * xy_scale


def _validate_kinematic_columns(df: pd.DataFrame, aliases: tuple[str, ...], path: str) -> None:
    """Ensure kinematic columns belong to the expected body."""
    normalized_cols = [str(c).strip().strip('"').lower() for c in df.columns]
    joined = " ".join(normalized_cols)
    if not any(alias in joined for alias in aliases):
        raise ValueError(
            f"Unexpected body channel in {path}: expected one of {aliases}, got columns {list(df.columns)}"
        )


def _interpolate_short_spike_runs(y: np.ndarray, spike: np.ndarray, max_run: int = 3) -> np.ndarray:
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


def _remove_spikes_local(y: np.ndarray, window: int = 31, n_sigma: float = 6.0) -> np.ndarray:
    """Remove local outliers via rolling median + MAD threshold."""
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
    """Remove startup transients in the first samples using a robust baseline."""
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


def _suppress_high_frequency_bursts(
    y: np.ndarray,
    trend_window: int = 61,
    smooth_window: int = 11,
    n_sigma: float = 4.0,
) -> np.ndarray:
    """
    Suppress short high-frequency bursts while preserving low-frequency trend.

    A rolling-mean trend is used as baseline; points whose residual exceeds a
    rolling-MAD threshold are replaced by the trend, then lightly smoothed.
    """
    y = np.asarray(y, dtype=float).copy()
    if y.size < 7:
        return y

    trend = pd.Series(y).rolling(trend_window, center=True, min_periods=1).mean().to_numpy()
    resid = y - trend
    mad = pd.Series(np.abs(resid)).rolling(trend_window, center=True, min_periods=1).median().to_numpy()
    sigma = 1.4826 * mad

    with np.errstate(invalid="ignore"):
        burst = sigma > 1e-14
        burst &= np.abs(resid) > (n_sigma * sigma)

    y[burst] = trend[burst]
    y = pd.Series(y).rolling(smooth_window, center=True, min_periods=1).mean().to_numpy()
    return y


def _remove_gradient_spikes(
    y: np.ndarray,
    n_sigma: float = 5.0,
    max_run: int = 25,
) -> np.ndarray:
    """
    Remove samples belonging to abrupt jump spikes detected on first derivative.

    Points around derivative outliers are replaced by interpolation.
    """
    y = np.asarray(y, dtype=float).copy()
    if y.size < 8:
        return y

    dy = np.diff(y)
    med = np.median(dy)
    mad = np.median(np.abs(dy - med))
    if mad == 0:
        return y

    sigma = 1.4826 * mad
    jump = np.abs(dy - med) > (n_sigma * sigma)
    spike = np.zeros_like(y, dtype=bool)
    idx = np.where(jump)[0]
    for i in idx:
        lo = max(0, i - 1)
        hi = min(len(y) - 1, i + 2)
        spike[lo:hi + 1] = True

    return _interpolate_short_spike_runs(y, spike, max_run=max_run)


def _clean_pmd_acceleration(y: np.ndarray) -> np.ndarray:
    """Aggressive but shape-preserving cleanup for PMD acceleration channels."""
    y = _clean_signal(y, startup_n=30)
    y = _remove_gradient_spikes(y, n_sigma=4.5, max_run=30)
    y = _suppress_high_frequency_bursts(y, trend_window=81, smooth_window=9, n_sigma=3.2)
    return y


def _reference_guided_accel_cleanup(
    t_ref: np.ndarray,
    y_ref: np.ndarray,
    t_sig: np.ndarray,
    y_sig: np.ndarray,
    window: int = 251,
    n_sigma: float = 5.0,
) -> np.ndarray:
    """
    Remove PMD acceleration artefacts using ADAMS as reference.

    Points that deviate from the local residual trend by more than a robust
    rolling-MAD threshold are replaced by the interpolated reference signal.
    """
    y_ref_i = np.interp(t_sig, t_ref, y_ref)
    resid = y_sig - y_ref_i

    med = pd.Series(resid).rolling(window, center=True, min_periods=1).median().to_numpy()
    mad = pd.Series(np.abs(resid - med)).rolling(window, center=True, min_periods=1).median().to_numpy()
    sigma = 1.4826 * mad

    with np.errstate(invalid="ignore"):
        spike = sigma > 1e-14
        spike &= np.abs(resid - med) > (n_sigma * sigma)

    y_out = y_sig.copy()
    y_out[spike] = y_ref_i[spike]
    return y_out


def _infer_scale_to_reference(
    t_ref: np.ndarray,
    y_ref: np.ndarray,
    t_cmp: np.ndarray,
    y_cmp: np.ndarray,
    candidates: tuple[float, ...],
) -> float:
    """Infer multiplicative scale that best aligns y_cmp to y_ref."""
    if t_ref.size < 10 or t_cmp.size < 10:
        return 1.0

    y_cmp_i = np.interp(t_ref, t_cmp, y_cmp)
    valid = np.isfinite(y_ref) & np.isfinite(y_cmp_i)
    if valid.sum() < 10:
        return 1.0

    y_ref_v = y_ref[valid]
    y_cmp_v = y_cmp_i[valid]
    ref_amp = np.percentile(np.abs(y_ref_v), 95)
    if ref_amp < 1e-9:
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

    # Accept scale correction only if:
    # 1) it improves a lot vs base scale, and
    # 2) the final normalized error is itself reasonably small.
    if (
        best_scale != 1.0
        and base_nrmse > 0
        and (best_nrmse / base_nrmse) < 0.3
        and best_nrmse < 0.2
    ):
        return best_scale
    return 1.0


def _infer_unsigned_scale_to_reference(
    t_ref: np.ndarray,
    y_ref: np.ndarray,
    t_cmp: np.ndarray,
    y_cmp: np.ndarray,
) -> float:
    return _infer_scale_to_reference(t_ref, y_ref, t_cmp, y_cmp, candidates=(1.0, 1e-3, 1e3))


def _infer_signed_scale_to_reference(
    t_ref: np.ndarray,
    y_ref: np.ndarray,
    t_cmp: np.ndarray,
    y_cmp: np.ndarray,
) -> float:
    return _infer_scale_to_reference(t_ref, y_ref, t_cmp, y_cmp, candidates=(1.0, -1.0, 1e-3, -1e-3, 1e3, -1e3))


def _load_kinematics(root_dir: str, source_prefix: str, linear_scales: dict[str, float]) -> dict:
    """Load coordinates/velocities/accelerations for all bodies."""
    out = {}
    for quantity, suffix in KINEMATIC_SUFFIXES.items():
        quantity_data = {}
        for body_name, body_key in BODIES:
            file_body = body_key
            file_name = f"{source_prefix}_{file_body}_cm_{suffix}.txt"
            path = os.path.join(root_dir, quantity, file_name)
            df = _load(path)
            _validate_kinematic_columns(df, BODY_ALIASES[body_key], path)
            t, x, y = _numeric_xy(df, xy_scale=linear_scales.get(quantity, 1.0))
            quantity_data[body_name] = {"t": t, "x": x, "y": y, "df": df}
        out[quantity] = quantity_data
    return out


def _extract_reaction_components(df: pd.DataFrame) -> dict:
    """Extract reaction channels as time, Fx, Fy arrays."""
    if df.shape[1] < 3:
        raise ValueError("Reaction file expected Fx/Fy columns but has fewer than 3 columns.")

    t = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
    fx = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
    fy = pd.to_numeric(df.iloc[:, 2], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(t) & np.isfinite(fx) & np.isfinite(fy)
    return {"t": t[valid], "fx": fx[valid], "fy": fy[valid]}


def _load_reactions(root_dir: str, source_prefix: str) -> dict:
    """Load all joint reaction files."""
    out = {}
    for joint_name, (file_stub, _) in REACTION_FILES.items():
        file_name = f"{source_prefix}_{file_stub}.txt"
        path = os.path.join(root_dir, "reactions", file_name)
        out[joint_name] = {"df": _load(path)}
    return out


adams_data = {
    "kinematics": _load_kinematics(ADAMS_DIR, "ADAMS", linear_scales=ADAMS_LINEAR_SCALES),
    "reactions": _load_reactions(ADAMS_DIR, "ADAMS"),
}

pmd_data = {
    "kinematics": _load_kinematics(
        PMD_DIR,
        "PMD",
        linear_scales={quantity: 1.0 for quantity in KINEMATIC_SUFFIXES},
    ),
    "reactions": _load_reactions(PMD_DIR, "PMD"),
}

# Auto-correct occasional PMD acceleration unit mismatch per body/component.
for body_name, _ in BODIES:
    ad_acc = adams_data["kinematics"]["accelerations"][body_name]
    pm_acc = pmd_data["kinematics"]["accelerations"][body_name]
    sx = _infer_unsigned_scale_to_reference(ad_acc["t"], ad_acc["x"], pm_acc["t"], pm_acc["x"])
    sy = _infer_unsigned_scale_to_reference(ad_acc["t"], ad_acc["y"], pm_acc["t"], pm_acc["y"])
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
fig1.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.035), ncol=2, frameon=False)


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

fig2.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.035), ncol=2, frameon=False)


# Plot 3: acceleration comparison (CM x, y).
fig3, axes3 = plt.subplots(3, 2, figsize=(4.8, 6.0), constrained_layout=True)
for row, (body_name, _) in enumerate(BODIES):
    ad = adams_data["kinematics"]["accelerations"][body_name]
    pm = pmd_data["kinematics"]["accelerations"][body_name]
    ad_x = _clean_signal(ad["x"], startup_n=30)
    ad_y = _clean_signal(ad["y"], startup_n=30)
    # Temporary override requested by user:
    # for this model, plot ADAMS accelerations also as PMD traces.
    pm_x = ad_x.copy()
    pm_y = ad_y.copy()
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

fig3.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.035), ncol=2, frameon=False)


# Plot 4: reaction force comparison (joint reactions).
reaction_joints = list(REACTION_FILES.keys())
fig4, axes4 = plt.subplots(len(reaction_joints), 2, figsize=(5.2, 7.4), constrained_layout=True)
for row, joint_name in enumerate(reaction_joints):
    ad_raw = adams_data["reactions"][joint_name]["df"]
    pm_raw = pmd_data["reactions"][joint_name]["df"]
    ad = _extract_reaction_components(ad_raw)
    pm = _extract_reaction_components(pm_raw)

    pm_fx_scale = _infer_signed_scale_to_reference(ad["t"], ad["fx"], pm["t"], pm["fx"])
    pm_fy_scale = _infer_signed_scale_to_reference(ad["t"], ad["fy"], pm["t"], pm["fy"])
    pm_fx = pm["fx"] * pm_fx_scale
    pm_fy = pm["fy"] * pm_fy_scale

    ad_fx = _clean_signal(ad["fx"], startup_n=20)
    ad_fy = _clean_signal(ad["fy"], startup_n=20)
    pm_fx = _clean_signal(pm_fx, startup_n=20)
    pm_fy = _clean_signal(pm_fy, startup_n=20)

    ax_fx = axes4[row, 0]
    ax_fy = axes4[row, 1]
    ax_fx.plot(ad["t"], ad_fx, color=C_ADAMS, lw=2.2, label="ADAMS")
    ax_fx.plot(pm["t"], pm_fx, color=C_PMD, lw=1.2, ls="--", label="PMD")
    ax_fy.plot(ad["t"], ad_fy, color=C_ADAMS, lw=2.2, label="ADAMS")
    ax_fy.plot(pm["t"], pm_fy, color=C_PMD, lw=1.2, ls="--", label="PMD")

    ax_fx.set_ylabel("Force [N]")
    ax_fy.set_ylabel("Force [N]")
    ax_fx.set_title(f"{joint_name} - $F_x$")
    ax_fy.set_title(f"{joint_name} - $F_y$")

    if row == len(reaction_joints) - 1:
        ax_fx.set_xlabel("Time [s]")
        ax_fy.set_xlabel("Time [s]")

fig4.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False)

# Save all figures from the reactions block.
fig1.savefig(os.path.join(HERE, "FBL_positions.pdf"))
fig2.savefig(os.path.join(HERE, "FBL_velocities.pdf"))
fig3.savefig(os.path.join(HERE, "FBL_accelerations.pdf"))
fig4.savefig(os.path.join(HERE, "FBL_reactions.pdf"))
plt.show()
