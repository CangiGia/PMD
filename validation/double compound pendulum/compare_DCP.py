"""
Double Compound Pendulum Validation — ADAMS vs pmd
===================================================
Compares CM positions and joint reaction forces produced by ADAMS/View
and pmd (CASADI-DAE) for the double compound pendulum.

Units
-----
  CM positions    : mm  (ADAMS) / m (pmd) → plotted in mm
  Forces          : N   (both ADAMS and pmd)
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ADAMS_DIR = os.path.join(HERE, 'ADAMS')
PMD_DIR   = os.path.join(HERE, 'PMD')

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
C_PMD   = '#2E7EC1'  # steel blue


def _load_pmd(path):
    return pd.read_csv(path, sep=r'\s+', header=0, engine='python')


def _load_adams(path):
    return pd.read_csv(path, sep=r'\s+', header=0, quotechar='"', engine='python')


A_cm = _load_adams(os.path.join(ADAMS_DIR, 'ADAMS_double_compound_pendulum_cm_coordinates.txt'))
A_rf = _load_adams(os.path.join(ADAMS_DIR, 'ADAMS_double_compound_pendulum_reaction_forces.txt'))

P_cm = _load_pmd(os.path.join(PMD_DIR, 'PMD_double_compound_pendulum_cm_coordinates.txt'))
P_rf = _load_pmd(os.path.join(PMD_DIR, 'PMD_double_compound_pendulum_reaction_forces.txt'))

# Figure 1 — CM Positions
links = ['Link 1', 'Link 2']

tA_cm = A_cm.iloc[:, 0].values
tP_cm = P_cm.iloc[:, 0].values

fig1, axes1 = plt.subplots(2, 2, figsize=(12, 6), constrained_layout=True)
fig1.suptitle('Double Compound Pendulum  —  CM Positions: ADAMS vs pmd',
              fontsize=13, fontweight='bold')

for i, link in enumerate(links):
    for j, coord in enumerate(['x', 'y']):
        ax = axes1[i, j]
        col = 2 * i + j + 1
        ax.plot(tA_cm, A_cm.iloc[:, col].values,          color=C_ADAMS, lw=4,   label='ADAMS')
        ax.plot(tP_cm, P_cm.iloc[:, col].values * 1e3,    color=C_PMD,   lw=1.5, ls='--', label='pmd')
        ax.set_ylabel(f'{coord}  [mm]')
        ax.set_title(f'{link}  —  {coord}')
        if i == 1:
            ax.set_xlabel('Time  [s]')
        ax.legend(loc='best')

# Figure 2 — Joint Reaction Forces
joints = ['Ground \u2013 Link 1', 'Link 1 \u2013 Link 2']

tA_rf = A_rf.iloc[:, 0].values
tP_rf = P_rf.iloc[:, 0].values

fig2, axes2 = plt.subplots(2, 2, figsize=(12, 6), constrained_layout=True)
fig2.suptitle('Double Compound Pendulum  —  Joint Reaction Forces: ADAMS vs pmd',
              fontsize=13, fontweight='bold')

for i, name in enumerate(joints):
    for j, coord in enumerate(['Fx', 'Fy']):
        ax = axes2[i, j]
        col = 2 * i + j + 1
        ax.plot(tA_rf, A_rf.iloc[:, col].values, color=C_ADAMS, lw=4,   label='ADAMS')
        ax.plot(tP_rf, P_rf.iloc[:, col].values, color=C_PMD,   lw=1.5, ls='--', label='pmd')
        ax.set_ylabel(f'{coord}  [N]')
        ax.set_title(f'{name}  —  {coord}')
        if i == 1:
            ax.set_xlabel('Time  [s]')
        ax.legend(loc='best')

plt.show()
