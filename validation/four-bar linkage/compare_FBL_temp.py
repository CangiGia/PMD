"""
Four-Bar Linkage Validation — ADAMS vs PMD  (presentation version)
===================================================================
Derived from compare_FBL.py with the following changes vs the original:

  1. Force spikes removed: isolated samples that deviate more than
     5× the inter-quartile range from the local median are replaced by
     linear interpolation, eliminating the vertical-stroke artefacts
     visible in some force subplots at t ≈ 0.
  2. No bold weight on any title or suptitle.
  3. Comic-Sans-like font (Comic Sans MS when available, fallback to
     cursive) applied to all text elements.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as _fm

HERE = os.path.dirname(os.path.abspath(__file__))
ADAMS_DIR = os.path.join(HERE, 'ADAMS')
PMD_DIR = os.path.join(HERE, 'PMD')

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


def _load_PMD(path):
    return pd.read_csv(path, sep=r'\s+', header=0, engine='python')


def _load_adams(path):
    return pd.read_csv(
        path, sep=r'\s+', header=0, skiprows=2, quotechar='"', engine='python'
    )


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
                lo = i - 1 if i > 0 else j
                hi = j     if j < n else i - 1
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
A_cm  = _load_adams(os.path.join(ADAMS_DIR, 'ADAMS_four_bar_linkage_cm_coordinates.txt'))
A_rf  = _load_adams(os.path.join(ADAMS_DIR, 'ADAMS_four_bar_linkage_reaction_forces.txt'))
A_mot = _load_adams(os.path.join(ADAMS_DIR, 'ADAMS_four_bar_linkage_motion.txt'))

P_cm  = _load_PMD(os.path.join(PMD_DIR, 'PMD_four_bar_linkage_cm_coordinates.txt'))
P_rf  = _load_PMD(os.path.join(PMD_DIR, 'PMD_four_bar_linkage_reaction_forces.txt'))
P_mot = _load_PMD(os.path.join(PMD_DIR, 'PMD_four_bar_linkage_motion.txt'))

# ------------------------------------------------------------------
# Figure 1 — CM Positions
# ------------------------------------------------------------------
links = ['Link 1 (Crank)', 'Link 2 (Coupler)', 'Link 3 (Follower)']

tA_cm = A_cm.iloc[:, 0].values
tP_cm = P_cm.iloc[:, 0].values

fig1, axes1 = plt.subplots(3, 2, figsize=(12, 9), constrained_layout=True)
fig1.suptitle('Four-Bar Linkage  —  CM Positions: ADAMS vs PMD',
              fontsize=13, fontweight='normal')

for i, link in enumerate(links):
    for j, coord in enumerate(['x', 'y']):
        ax = axes1[i, j]
        col = 2 * i + j + 1  # columns: 1=link1_x, 2=link1_y, 3=link2_x, ...
        ax.plot(tA_cm, A_cm.iloc[:, col].values, color=C_ADAMS, lw=4, label='ADAMS')
        ax.plot(tP_cm, P_cm.iloc[:, col].values, color=C_PMD, lw=1.5, ls='--', label='PMD')
        ax.set_ylabel(f'{coord}  [mm]')
        ax.set_title(f'{link}  —  {coord}', fontweight='normal')
        if i == 2:
            ax.set_xlabel('Time  [s]')
        ax.legend(loc='best')

# ------------------------------------------------------------------
# Figure 2 — Joint Reaction Forces
# ------------------------------------------------------------------
joints = [
    'Ground \u2013 Link 1',
    'Link 1 \u2013 Link 2',
    'Link 2 \u2013 Link 3',
    'Link 3 \u2013 Ground',
]

tA_rf = A_rf.iloc[:, 0].values
tP_rf = P_rf.iloc[:, 0].values

fig2, axes2 = plt.subplots(4, 2, figsize=(12, 12), constrained_layout=True)
fig2.suptitle('Four-Bar Linkage  —  Joint Reaction Forces: ADAMS vs PMD',
              fontsize=13, fontweight='normal')

# Fix: skip the first sample (t=0 initial-condition transient) and remove
# isolated spikes to eliminate vertical-stroke artefacts.
S = 1  # number of leading samples to drop

for i, name in enumerate(joints):
    for j, (coord, lbl) in enumerate([('Fx', 'Fx  [N]'), ('Fy', 'Fy  [N]')]):
        ax = axes2[i, j]
        col = 2 * i + j + 1
        fA = _remove_spikes(A_rf.iloc[S:, col].values)
        fP = _remove_spikes(P_rf.iloc[S:, col].values)
        ax.plot(tA_rf[S:], fA, color=C_ADAMS, lw=4, label='ADAMS')
        ax.plot(tP_rf[S:], fP, color=C_PMD, lw=1.5, ls='--', label='PMD')
        ax.set_ylabel(lbl)
        ax.set_title(f'{name}  —  {coord}', fontweight='normal')
        if i == 3:
            ax.set_xlabel('Time  [s]')
        ax.legend(loc='best')

# ------------------------------------------------------------------
# Figure 3 — Motion Driver Torque
#   ADAMS: Element_Torque.Z in N·mm  →  convert to N·m (÷ 1000)
#   PMD:   lambda_0 (Lagrange multiplier for motion constraint) in N·m
# ------------------------------------------------------------------
tA_mot = A_mot.iloc[:, 0].values
torq_A = _remove_spikes(A_mot.iloc[S:, 1].values) * 1e-3  # N·mm → N·m

tP_mot = P_mot.iloc[:, 0].values
torq_P = _remove_spikes(P_mot.iloc[S:, 1].values)  # lambda_0 [N·m]

fig3, ax3 = plt.subplots(figsize=(9, 4), constrained_layout=True)
fig3.suptitle('Four-Bar Linkage  —  Motion Driver Torque: ADAMS vs PMD',
              fontsize=13, fontweight='normal')

ax3.plot(tA_mot[S:], torq_A, color=C_ADAMS, lw=4, label='ADAMS  (torque Z)')
ax3.plot(tP_mot[S:], torq_P, color=C_PMD, lw=1.5, ls='--',
         label=r'PMD  ($\lambda_0$, motion constraint)')
ax3.set_ylabel('Torque  [N\u00b7m]')
ax3.set_xlabel('Time  [s]')
ax3.legend(loc='best')

fig1.savefig(os.path.join(HERE, 'FBL_positions.png'),    dpi=500, bbox_inches='tight')
fig2.savefig(os.path.join(HERE, 'FBL_forces.png'),       dpi=500, bbox_inches='tight')

plt.show()
