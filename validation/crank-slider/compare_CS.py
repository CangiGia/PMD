"""
Crank-Slider Validation — ADAMS vs PMD
=======================================
Compares CM positions and joint reaction forces produced by ADAMS/View
and PMD (CASADI-DAE) for the crank-slider mechanism.

Units
-----
  ADAMS CM positions : mm  →  converted to m for plotting
  PMD   CM positions : m
  Forces (both)      : N
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ADAMS_DIR = os.path.join(HERE, 'ADAMS')
PMD_DIR = os.path.join(HERE, 'PMD')

# ── Style ───────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'figure.facecolor': '#ffffff',
    'axes.facecolor': '#f5f5f5',
    'axes.edgecolor': '#cccccc',
    'axes.grid': True,
    'grid.color': '#e0e0e0',
    'grid.linewidth': 0.8,
    'grid.linestyle': '--',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'axes.titleweight': 'semibold',
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'legend.frameon': True,
    'legend.framealpha': 0.85,
    'legend.edgecolor': '#cccccc',
    'legend.fontsize': 8,
    'lines.linewidth': 1.8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
})

C_ADAMS = '#D64B35'  # warm red
C_PMD = '#2E7EC1'  # steel blue


def _load(path):
    """Load a space/tab-separated validation file with one quoted header row."""
    return pd.read_csv(path, sep=r'\s+', header=0, quotechar='"', engine='python')


# ── CM positions ─────────────────────────────────────────────────────────────
A_crank_cm = _load(os.path.join(ADAMS_DIR, 'ADAMS_crankshaft_cm_coordinates.txt'))
A_rod_cm = _load(os.path.join(ADAMS_DIR, 'ADAMS_rod_cm_coordinates.txt'))
A_slider_cm = _load(os.path.join(ADAMS_DIR, 'ADAMS_slider_cm_coordinates.txt'))

P_crank_cm = _load(os.path.join(PMD_DIR, 'PMD_crankshaft_cm_coordinates.txt'))
P_rod_cm = _load(os.path.join(PMD_DIR, 'PMD_rod_cm_coordinates.txt'))
P_slider_cm = _load(os.path.join(PMD_DIR, 'PMD_slider_cm_coordinates.txt'))

# ADAMS: mm → m
tA_cm = [d.iloc[:, 0].values for d in (A_crank_cm, A_rod_cm, A_slider_cm)]
xA = [d.iloc[:, 1].values * 1e-3 for d in (A_crank_cm, A_rod_cm, A_slider_cm)]
yA = [d.iloc[:, 2].values * 1e-3 for d in (A_crank_cm, A_rod_cm, A_slider_cm)]

tP_cm = [d.iloc[:, 0].values for d in (P_crank_cm, P_rod_cm, P_slider_cm)]
xP = [d.iloc[:, 1].values for d in (P_crank_cm, P_rod_cm, P_slider_cm)]
yP = [d.iloc[:, 2].values for d in (P_crank_cm, P_rod_cm, P_slider_cm)]

# ── Joint reaction forces ─────────────────────────────────────────────────
A_crank_gnd = _load(os.path.join(ADAMS_DIR, 'ADAMS_crankshaft_ground_forces.txt'))
A_crank_rod = _load(os.path.join(ADAMS_DIR, 'ADAMS_crankshaft_rod_forces.txt'))
A_rod_sldr = _load(os.path.join(ADAMS_DIR, 'ADAMS_rod_slider_forces.txt'))
A_sldr_gnd = _load(os.path.join(ADAMS_DIR, 'ADAMS_slider_ground_forces.txt'))

P_crank_gnd = _load(os.path.join(PMD_DIR, 'PMD_crankshaft_ground_forces.txt'))
P_crank_rod = _load(os.path.join(PMD_DIR, 'PMD_crankshaft_rod_forces.txt'))
P_rod_sldr = _load(os.path.join(PMD_DIR, 'PMD_rod_slider_forces.txt'))
P_sldr_gnd = _load(os.path.join(PMD_DIR, 'PMD_slider_ground_forces.txt'))

# ── Figure 1 — CM Positions ──────────────────────────────────────────────
bodies = ['Crankshaft', 'Rod', 'Slider']

fig1, axes1 = plt.subplots(3, 2, figsize=(12, 9), constrained_layout=True)
fig1.suptitle('Crank-Slider  —  CM Positions: ADAMS vs PMD',
              fontsize=13, fontweight='bold')

for i, body in enumerate(bodies):
    for j, (coord, col_A, col_P) in enumerate([
        ('x  [m]', xA[i], xP[i]),
        ('y  [m]', yA[i], yP[i]),
    ]):
        ax = axes1[i, j]
        ax.plot(tA_cm[i], col_A, color=C_ADAMS, lw=1.8, label='ADAMS')
        ax.plot(tP_cm[i], col_P, color=C_PMD, lw=1.4, ls='--', label='PMD')
        ax.set_ylabel(coord)
        ax.set_title(f'{body}  —  {"x" if j == 0 else "y"}')
        if i == 2:
            ax.set_xlabel('Time  [s]')
        ax.legend(loc='best')

# ── Figure 2 — Joint Reaction Forces ─────────────────────────────────────
#   Slider–Ground: ADAMS exports (Fx≈0, Fy); PMD exports only F_perp (= Fy).
#   The left subplot is hidden for that row.
joints = [
    ('Ground – Crankshaft', A_crank_gnd, P_crank_gnd, True),
    ('Crankshaft – Rod', A_crank_rod, P_crank_rod, True),
    ('Rod – Slider', A_rod_sldr, P_rod_sldr, True),
    ('Slider – Ground', A_sldr_gnd, P_sldr_gnd, False),
]

fig2, axes2 = plt.subplots(4, 2, figsize=(12, 12), constrained_layout=True)
fig2.suptitle('Crank-Slider  —  Joint Reaction Forces: ADAMS vs PMD',
              fontsize=13, fontweight='bold')

for i, (name, A_df, P_df, has_fx) in enumerate(joints):
    ax_l = axes2[i, 0]
    ax_r = axes2[i, 1]
    tA = A_df.iloc[:, 0].values
    tP = P_df.iloc[:, 0].values

    if has_fx:
        ax_l.plot(tA, A_df.iloc[:, 1].values, color=C_ADAMS, lw=1.8, label='ADAMS')
        ax_l.plot(tP, P_df.iloc[:, 1].values, color=C_PMD, lw=1.4, ls='--', label='PMD')
        ax_l.set_ylabel('Fx  [N]')
        ax_l.set_title(f'{name}  —  Fx')
        ax_l.legend(loc='best')

        ax_r.plot(tA, A_df.iloc[:, 2].values, color=C_ADAMS, lw=1.8, label='ADAMS')
        ax_r.plot(tP, P_df.iloc[:, 2].values, color=C_PMD, lw=1.4, ls='--', label='PMD')
        ax_r.set_ylabel('Fy  [N]')
        ax_r.set_title(f'{name}  —  Fy')
        ax_r.legend(loc='best')
    else:
        # PMD exports only F_perp (normal to sliding direction = Fy)
        ax_l.set_visible(False)
        ax_r.plot(tA, A_df.iloc[:, 2].values, color=C_ADAMS, lw=1.8, label='ADAMS  Fy')
        ax_r.plot(tP, P_df.iloc[:, 1].values, color=C_PMD, lw=1.4, ls='--', label='PMD  F\u22a5')
        ax_r.set_ylabel('F  [N]')
        ax_r.set_title(f'{name}  —  F\u22a5 / Fy')
        ax_r.legend(loc='best')

    if i == 3:
        ax_r.set_xlabel('Time  [s]')

plt.show()
