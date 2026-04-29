"""
Mass-Spring-Damper Validation — ADAMS vs PMD
=============================================
Compares CM displacement of the mass produced by ADAMS/View and PMD
(CASADI-DAE) for the mass-spring-damper system.

Model parameters
----------------
  mass       : 1000 kg
  stiffness  : 1 000 000 N/m
  damping    : 6 500 N·s/m
  v0         : 10 m/s  (upward)
  L0         : 1.5 m   (no preload at t=0)
  gravity    : absent

Units
-----
  CM displacement (Y) : m   (both ADAMS and PMD)
  Time                : s   (both ADAMS and PMD)
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ADAMS_DIR = os.path.join(HERE, 'ADAMS')
PMD_DIR = os.path.join(HERE, 'PMD')

# ── Style ────────────────────────────────────────────────────────────────
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
    return pd.read_csv(path, sep='\t', header=0, engine='python')


def _load_adams(path):
    return pd.read_csv(path, sep=r'\s+', header=0, quotechar='"', engine='python')


# ── Load data ─────────────────────────────────────────────────────────────
A_cm = _load_adams(os.path.join(ADAMS_DIR, 'ADAMS_mass_spring_damper_cm_coordinates.txt'))
P_cm = _load_pmd(os.path.join(PMD_DIR, 'PMD_mass_spring_damper_cm_coordinates.txt'))

# Column aliases
tA = A_cm.iloc[:, 0].values          # Time [s]
yA = A_cm.iloc[:, 2].values          # CM_Position.Y [m]

tP = P_cm.iloc[:, 0].values           # time [s]
yP = P_cm.iloc[:, 2].values           # CM displacement Y [m]  (col: 'mass / y')  (col index 2)

# ── Figure — CM Displacement Y ────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
fig.suptitle('Mass-Spring-Damper  —  CM Displacement Y: ADAMS vs PMD',
             fontsize=13, fontweight='bold')

ax.plot(tA, yA, color=C_ADAMS, lw=4,   label='ADAMS')
ax.plot(tP, yP, color=C_PMD,   lw=1.5, ls='--', label='PMD')
ax.set_xlabel('Time  [s]')
ax.set_ylabel('CM displacement Y  [m]')
ax.set_title('mass  —  y')
ax.legend(loc='best')

plt.show()
