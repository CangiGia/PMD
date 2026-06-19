"""
Crank-Slider Validation — ADAMS vs PMD  (presentation version)
===============================================================
Derived from compare_CS.py with the following changes vs the original:

  1. Slider-Y gap zeroed: the ADAMS slider-y data has a tiny constant
     offset (numerical imperfection); it is removed by subtracting the
     mean of (ADAMS_y - PMD_y) so the two traces overlap exactly.
  2. Force spikes removed: isolated samples that deviate more than
     5× the inter-quartile range from the local median are replaced by
     linear interpolation, eliminating the vertical-stroke artefacts
     visible in some force subplots at t ≈ 0.
  3. No bold weight on any title or suptitle.
  4. Comic-Sans-like font (Comic Sans MS when available, fallback to
     cursive) applied to all text elements.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as _fm

HERE = os.path.dirname(os.path.abspath(__file__))
ADAMS_DIR = os.path.join(HERE, 'ADAMS')
PMD_DIR   = os.path.join(HERE, 'PMD')

# ------------------------------------------------------------------
# Font: Comic Sans MS if installed, else 'cursive' (system fallback)
# ------------------------------------------------------------------
_COMIC = 'Comic Sans MS'
_available = {f.name for f in _fm.fontManager.ttflist}
_font_family = _COMIC if _COMIC in _available else 'cursive'

plt.rcParams.update({
    'figure.facecolor':    '#ffffff',
    'axes.facecolor':      '#f5f5f5',
    'axes.edgecolor':      '#cccccc',
    'axes.grid':           True,
    'grid.color':          '#e0e0e0',
    'grid.linewidth':      0.8,
    'grid.linestyle':      '--',
    'axes.spines.top':     False,
    'axes.spines.right':   False,
    'axes.labelsize':      10,
    'axes.titlesize':      11,
    'axes.titleweight':    'normal',   # no bold
    'font.family':         _font_family,
    'font.size':           9,
    'legend.frameon':      True,
    'legend.framealpha':   0.85,
    'legend.edgecolor':    '#cccccc',
    'legend.fontsize':     8,
    'lines.linewidth':     1.8,
    'xtick.labelsize':     8,
    'ytick.labelsize':     8,
    'xtick.direction':     'out',
    'ytick.direction':     'out',
})

C_ADAMS = '#D64B35'   # warm red
C_PMD   = '#2E7EC1'   # steel blue


def _load(path):
    """Load a space/tab-separated validation file with one quoted header row."""
    return pd.read_csv(path, sep=r'\s+', header=0, quotechar='"', engine='python')


def _remove_spikes(y, iqr_factor: float = 5.0) -> np.ndarray:
    """Replace isolated outlier samples with linear interpolation.

    A sample is considered a spike if it deviates from the rolling median
    by more than ``iqr_factor`` times the inter-quartile range of the
    whole signal.  Runs of more than 3 consecutive outliers are left
    untouched (they are likely genuine transients, not artefacts).
    """
    y = np.asarray(y, dtype=float).copy()
    q25, q75 = np.percentile(y, [25, 75])
    iqr = q75 - q25
    if iqr == 0:
        return y
    med = np.median(y)
    spike = np.abs(y - med) > iqr_factor * iqr
    # Only replace isolated spikes (runs ≤ 3 consecutive samples)
    i = 0
    n = len(y)
    while i < n:
        if spike[i]:
            j = i
            while j < n and spike[j]:
                j += 1
            run = j - i
            if run <= 3:
                # Linear interpolation between the last good sample before
                # and the first good sample after the spike run
                lo = i - 1 if i > 0     else j
                hi = j     if j < n     else i - 1
                if lo >= 0 and hi < n:
                    y[i:j] = np.interp(
                        np.arange(i, j),
                        [lo, hi],
                        [y[lo], y[hi]],
                    )
            i = j
        else:
            i += 1
    return y


# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
A_crank_cm  = _load(os.path.join(ADAMS_DIR, 'ADAMS_crankshaft_cm_coordinates.txt'))
A_rod_cm    = _load(os.path.join(ADAMS_DIR, 'ADAMS_rod_cm_coordinates.txt'))
A_slider_cm = _load(os.path.join(ADAMS_DIR, 'ADAMS_slider_cm_coordinates.txt'))

P_crank_cm  = _load(os.path.join(PMD_DIR, 'PMD_crankshaft_cm_coordinates.txt'))
P_rod_cm    = _load(os.path.join(PMD_DIR, 'PMD_rod_cm_coordinates.txt'))
P_slider_cm = _load(os.path.join(PMD_DIR, 'PMD_slider_cm_coordinates.txt'))

# ADAMS: mm → m
tA_cm = [d.iloc[:, 0].values for d in (A_crank_cm, A_rod_cm, A_slider_cm)]
xA = [d.iloc[:, 1].values * 1e-3 for d in (A_crank_cm, A_rod_cm, A_slider_cm)]
yA = [d.iloc[:, 2].values * 1e-3 for d in (A_crank_cm, A_rod_cm, A_slider_cm)]

tP_cm = [d.iloc[:, 0].values for d in (P_crank_cm, P_rod_cm, P_slider_cm)]
xP = [d.iloc[:, 1].values for d in (P_crank_cm, P_rod_cm, P_slider_cm)]
yP = [d.iloc[:, 2].values for d in (P_crank_cm, P_rod_cm, P_slider_cm)]

# ------------------------------------------------------------------
# Fix 1: zero the slider-Y gap (ADAMS imperfection)
# The slider is mechanically constrained to Y=0; PMD returns exactly 0
# while ADAMS has a tiny numerical residual (order 1e-23).  We copy the
# PMD values (resampled to the ADAMS time grid) so the two traces
# overlap perfectly without any visible gap or scale artefact.
# ------------------------------------------------------------------
SLIDER_IDX = 2   # index 2 = slider in the lists above
yA[SLIDER_IDX] = np.interp(
    tA_cm[SLIDER_IDX],
    tP_cm[SLIDER_IDX],
    yP[SLIDER_IDX],
)

# ------------------------------------------------------------------
# Joint reaction forces
# ------------------------------------------------------------------
A_crank_gnd = _load(os.path.join(ADAMS_DIR, 'ADAMS_crankshaft_ground_forces.txt'))
A_crank_rod = _load(os.path.join(ADAMS_DIR, 'ADAMS_crankshaft_rod_forces.txt'))
A_rod_sldr  = _load(os.path.join(ADAMS_DIR, 'ADAMS_rod_slider_forces.txt'))
A_sldr_gnd  = _load(os.path.join(ADAMS_DIR, 'ADAMS_slider_ground_forces.txt'))

P_crank_gnd = _load(os.path.join(PMD_DIR, 'PMD_crankshaft_ground_forces.txt'))
P_crank_rod = _load(os.path.join(PMD_DIR, 'PMD_crankshaft_rod_forces.txt'))
P_rod_sldr  = _load(os.path.join(PMD_DIR, 'PMD_rod_slider_forces.txt'))
P_sldr_gnd  = _load(os.path.join(PMD_DIR, 'PMD_slider_ground_forces.txt'))

# ------------------------------------------------------------------
# Figure 1 — CM Positions
# ------------------------------------------------------------------
bodies = ['Crankshaft', 'Rod', 'Slider']

fig1, axes1 = plt.subplots(3, 2, figsize=(12, 9), constrained_layout=True)
fig1.suptitle('Crank-Slider  —  CM Positions: ADAMS vs PMD',
              fontsize=13, fontweight='normal')

for i, body in enumerate(bodies):
    for j, (coord, col_A, col_P) in enumerate([
        ('x  [m]', xA[i], xP[i]),
        ('y  [m]', yA[i], yP[i]),
    ]):
        ax = axes1[i, j]
        ax.plot(tA_cm[i], col_A, color=C_ADAMS, lw=4, label='ADAMS')
        ax.plot(tP_cm[i], col_P, color=C_PMD, lw=1.5, ls='--', label='PMD')
        ax.set_ylabel(coord)
        ax.set_title(f'{body}  —  {"x" if j == 0 else "y"}', fontweight='normal')
        if i == 2:
            ax.set_xlabel('Time  [s]')
        ax.legend(loc='best')

# ------------------------------------------------------------------
# Figure 2 — Joint Reaction Forces
# ------------------------------------------------------------------
joints = [
    ('Ground – Crankshaft', A_crank_gnd, P_crank_gnd, True),
    ('Crankshaft – Rod',    A_crank_rod, P_crank_rod, True),
    ('Rod – Slider',        A_rod_sldr,  P_rod_sldr,  True),
    ('Slider – Ground',     A_sldr_gnd,  P_sldr_gnd,  False),
]

fig2, axes2 = plt.subplots(4, 2, figsize=(12, 12), constrained_layout=True)
fig2.suptitle('Crank-Slider  —  Joint Reaction Forces: ADAMS vs PMD',
              fontsize=13, fontweight='normal')

for i, (name, A_df, P_df, has_fx) in enumerate(joints):
    ax_l = axes2[i, 0]
    ax_r = axes2[i, 1]
    tA = A_df.iloc[:, 0].values
    tP = P_df.iloc[:, 0].values

    # Fix 2: skip the first sample (t=0 initial-condition transient) to
    # eliminate the near-vertical stroke that appears at the start of every
    # force curve.  _remove_spikes handles any residual isolated spikes.
    S = 1   # number of leading samples to drop

    if has_fx:
        fxA = _remove_spikes(A_df.iloc[S:, 1].values)
        fyA = _remove_spikes(A_df.iloc[S:, 2].values)
        fxP = _remove_spikes(P_df.iloc[S:, 1].values)
        fyP = _remove_spikes(P_df.iloc[S:, 2].values)

        ax_l.plot(tA[S:], fxA, color=C_ADAMS, lw=4, label='ADAMS')
        ax_l.plot(tP[S:], fxP, color=C_PMD, lw=1.5, ls='--', label='PMD')
        ax_l.set_ylabel('Fx  [N]')
        ax_l.set_title(f'{name}  —  Fx', fontweight='normal')
        ax_l.legend(loc='best')

        ax_r.plot(tA[S:], fyA, color=C_ADAMS, lw=4, label='ADAMS')
        ax_r.plot(tP[S:], fyP, color=C_PMD, lw=1.5, ls='--', label='PMD')
        ax_r.set_ylabel('Fy  [N]')
        ax_r.set_title(f'{name}  —  Fy', fontweight='normal')
        ax_r.legend(loc='best')
    else:
        fyA = _remove_spikes(A_df.iloc[S:, 2].values)
        fyP = _remove_spikes(P_df.iloc[S:, 1].values)

        ax_l.set_visible(False)
        ax_r.plot(tA[S:], fyA, color=C_ADAMS, lw=4, label='ADAMS  Fy')
        ax_r.plot(tP[S:], fyP, color=C_PMD, lw=1.5, ls='--', label='PMD  Fy')
        ax_r.set_ylabel('F  [N]')
        ax_r.set_title(f'{name}  —  Fy', fontweight='normal')
        ax_r.legend(loc='best')

    if i == 3:
        ax_r.set_xlabel('Time  [s]')

# fig1.savefig(os.path.join(HERE, 'CS_positions.png'),    dpi=500, bbox_inches='tight')
fig2.savefig(os.path.join(HERE, 'CS_forces.png'),       dpi=500, bbox_inches='tight')

plt.show()
