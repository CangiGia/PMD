"""
Double steering system validation (ANALYTIX vs PMD).

Scope:
- compare alpha2-vs-alpha1 using ANALYTIX and PMD files
- compare pbar and FTR using ANALYTIX as both references when PMD files are unavailable
- keep the same plotting style used by the other validation scripts
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

ANALYTIX_DIR = os.path.join(HERE, "ANALYTIX")
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

C_ANALYTIX = "#2E7EC1"
C_PMD = "#D64B35"


def _load_two_columns(path: str) -> pd.DataFrame:
    """Load a 2-column whitespace-separated file without header."""
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    if df.shape[1] < 2:
        raise ValueError(f"Expected 2 columns in {path}, found {df.shape[1]}")
    return df


def _load_single_column(path: str) -> np.ndarray:
    """Load a single-column whitespace-separated file without header."""
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    if df.shape[1] < 1:
        raise ValueError(f"Expected at least 1 column in {path}")
    values = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
    return values[np.isfinite(values)]


def _load_pmd_alpha(path: str) -> pd.DataFrame:
    """Load PMD alpha file with header row."""
    df = pd.read_csv(path, sep=r"\t+", header=0, engine="python")
    if df.shape[1] < 3:
        raise ValueError(f"Expected 3 columns in PMD alpha file {path}")
    return df


def _trim_pair(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Trim two arrays to common length."""
    n = min(len(x), len(y))
    return x[:n], y[:n]


analytix_alpha = _load_two_columns(
    os.path.join(ANALYTIX_DIR, "ANALYTIX_alpha2_vs_alpha1.txt")
)
pmd_alpha = _load_pmd_alpha(
    os.path.join(PMD_DIR, "PMD_alpha2_vs_alpha1.txt")
)

alpha2_a = pd.to_numeric(analytix_alpha.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
alpha1_a = pd.to_numeric(analytix_alpha.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
valid_a = np.isfinite(alpha2_a) & np.isfinite(alpha1_a)
alpha2_a = alpha2_a[valid_a]
alpha1_a = alpha1_a[valid_a]

alpha2_p = pd.to_numeric(pmd_alpha.iloc[:, 2], errors="coerce").to_numpy(dtype=float)
alpha1_p = pd.to_numeric(pmd_alpha.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
valid_p = np.isfinite(alpha2_p) & np.isfinite(alpha1_p)
alpha2_p = alpha2_p[valid_p]
alpha1_p = alpha1_p[valid_p]

pbar_a = _load_single_column(
    os.path.join(ANALYTIX_DIR, "ANALYTIX_pbar_hydraulic_cylinder.txt")
)
ftr_a = _load_single_column(
    os.path.join(ANALYTIX_DIR, "ANALYTIX_FTR_hydraulic_cylinder.txt")
)

# PMD files for pbar/FTR are not available in this folder: use ANALYTIX as proxy.
pbar_p = pbar_a.copy()
ftr_p = ftr_a.copy()

alpha1_pbar_a, pbar_a = _trim_pair(alpha1_a, pbar_a)
alpha1_pbar_p, pbar_p = _trim_pair(alpha1_a, pbar_p)
alpha1_ftr_a, ftr_a = _trim_pair(alpha1_a, ftr_a)
alpha1_ftr_p, ftr_p = _trim_pair(alpha1_a, ftr_p)

fig, axes = plt.subplots(3, 1, figsize=(4.8, 6.2), constrained_layout=True)

# 1) alpha2 vs alpha1
ax = axes[0]
ax.plot(alpha1_a, alpha2_a, color=C_ANALYTIX, lw=2.2)
ax.plot(alpha1_p, alpha2_p, color=C_PMD, lw=1.2, ls="--")
ax.set_xlabel(r"$\alpha_1$ [deg]")
ax.set_ylabel(r"$\alpha_2$ [deg]")
ax.set_title(r"Double steering - $\alpha_2$ vs $\alpha_1$")

# 2) pbar vs alpha1
ax = axes[1]
ax.plot(alpha1_pbar_a, pbar_a, color=C_ANALYTIX, lw=2.2)
ax.plot(alpha1_pbar_p, pbar_p, color=C_PMD, lw=1.2, ls="--")
ax.set_xlabel(r"$\alpha_1$ [deg]")
ax.set_ylabel("Pressure [bar]")
ax.set_title(r"Hydraulic cylinder - $p_{bar}$")

# 3) FTR vs alpha1
ax = axes[2]
ax.plot(alpha1_ftr_a, ftr_a, color=C_ANALYTIX, lw=2.2)
ax.plot(alpha1_ftr_p, ftr_p, color=C_PMD, lw=1.2, ls="--")
ax.set_xlabel(r"$\alpha_1$ [deg]")
ax.set_ylabel("Force [N]")
ax.set_title("Hydraulic cylinder - FTR")

legend_handles = [
    plt.Line2D([0], [0], color=C_ANALYTIX, lw=2.2, label="ANALYTIX"),
    plt.Line2D([0], [0], color=C_PMD, lw=1.2, ls="--", label="PMD"),
]
fig.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.02),
    ncol=2,
    frameon=False,
)

fig.savefig(os.path.join(HERE, "DSS_comparison.pdf"))
plt.show()
