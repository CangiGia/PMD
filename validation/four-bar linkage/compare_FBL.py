"""
Four-Bar Linkage Validation — ADAMS vs pmd
==========================================
Compares CM positions, joint reaction forces, and motion driver torque
produced by ADAMS/View and pmd (CASADI-DAE) for the four-bar linkage.

Units
-----
  CM positions    : mm  (both ADAMS and pmd)
  Forces          : N   (both ADAMS and pmd)
  Motion torque   : ADAMS N·mm  →  converted to N·m to match pmd (lambda_0)

ADAMS file format: 2 extra header lines (model title + blank) before column names.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ADAMS_DIR = os.path.join(HERE, 'ADAMS')
PMD_DIR = os.path.join(HERE, 'pmd')

# Style
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


def _load_pmd(path):
    return pd.read_csv(path, sep=r'\s+', header=0, engine='python')


def _load_adams(path):
    return pd.read_csv(
        path, sep=r'\s+', header=0, skiprows=2, quotechar='"', engine='python'
    )


# Load data
A_cm = _load_adams(os.path.join(ADAMS_DIR, 'ADAMS_four_bar_linkage_cm_coordinates.txt'))
A_rf = _load_adams(os.path.join(ADAMS_DIR, 'ADAMS_four_bar_linkage_reaction_forces.txt'))
A_mot = _load_adams(os.path.join(ADAMS_DIR, 'ADAMS_four_bar_linkage_motion.txt'))

P_cm = _load_pmd(os.path.join(PMD_DIR, 'PMD_four_bar_linkage_cm_coordinates.txt'))
P_rf = _load_pmd(os.path.join(PMD_DIR, 'PMD_four_bar_linkage_reaction_forces.txt'))
P_mot = _load_pmd(os.path.join(PMD_DIR, 'PMD_four_bar_linkage_motion.txt'))

# Figure 1 — CM Positions
links = ['Link 1 (Crank)', 'Link 2 (Coupler)', 'Link 3 (Follower)']

tA_cm = A_cm.iloc[:, 0].values
tP_cm = P_cm.iloc[:, 0].values

fig1, axes1 = plt.subplots(3, 2, figsize=(12, 9), constrained_layout=True)
fig1.suptitle('Four-Bar Linkage  —  CM Positions: ADAMS vs pmd',
              fontsize=13, fontweight='bold')

for i, link in enumerate(links):
    for j, coord in enumerate(['x', 'y']):
        ax = axes1[i, j]
        col = 2 * i + j + 1  # columns: 1=link1_x, 2=link1_y, 3=link2_x, ...
        ax.plot(tA_cm, A_cm.iloc[:, col].values, color=C_ADAMS, lw=4, label='ADAMS')
        ax.plot(tP_cm, P_cm.iloc[:, col].values, color=C_PMD, lw=1.5, ls='--', label='pmd')
        ax.set_ylabel(f'{coord}  [mm]')
        ax.set_title(f'{link}  —  {coord}')
        if i == 2:
            ax.set_xlabel('Time  [s]')
        ax.legend(loc='best')

# Figure 2 — Joint Reaction Forces
joints = [
    'Ground \u2013 Link 1',
    'Link 1 \u2013 Link 2',
    'Link 2 \u2013 Link 3',
    'Link 3 \u2013 Ground',
]

tA_rf = A_rf.iloc[:, 0].values
tP_rf = P_rf.iloc[:, 0].values

fig2, axes2 = plt.subplots(4, 2, figsize=(12, 12), constrained_layout=True)
fig2.suptitle('Four-Bar Linkage  —  Joint Reaction Forces: ADAMS vs pmd',
              fontsize=13, fontweight='bold')

for i, name in enumerate(joints):
    for j, (coord, lbl) in enumerate([('Fx', 'Fx  [N]'), ('Fy', 'Fy  [N]')]):
        ax = axes2[i, j]
        col = 2 * i + j + 1
        ax.plot(tA_rf, A_rf.iloc[:, col].values, color=C_ADAMS, lw=4, label='ADAMS')
        ax.plot(tP_rf, P_rf.iloc[:, col].values, color=C_PMD, lw=1.5, ls='--', label='pmd')
        ax.set_ylabel(lbl)
        ax.set_title(f'{name}  —  {coord}')
        if i == 3:
            ax.set_xlabel('Time  [s]')
        ax.legend(loc='best')

# Figure 3 — Motion Driver Torque
#   ADAMS: Element_Torque.Z in N·mm  →  convert to N·m (÷ 1000)
#   pmd:   lambda_0 (Lagrange multiplier for motion constraint) in N·m
tA_mot = A_mot.iloc[:, 0].values
torq_A = A_mot.iloc[:, 1].values * 1e-3  # N·mm → N·m

tP_mot = P_mot.iloc[:, 0].values
torq_P = P_mot.iloc[:, 1].values  # lambda_0 [N·m]

fig3, ax3 = plt.subplots(figsize=(9, 4), constrained_layout=True)
fig3.suptitle('Four-Bar Linkage  —  Motion Driver Torque: ADAMS vs pmd',
              fontsize=13, fontweight='bold')

ax3.plot(tA_mot, torq_A, color=C_ADAMS, lw=4, label='ADAMS  (torque Z)')
ax3.plot(tP_mot, torq_P, color=C_PMD, lw=1.5, ls='--',
         label=r'pmd  ($\lambda_0$, motion constraint)')
ax3.set_ylabel('Torque  [N\u00b7m]')
ax3.set_xlabel('Time  [s]')
ax3.legend(loc='best')

plt.show()
