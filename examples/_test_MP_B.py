"""
MP_B Model - Planar Multibody Dynamics
========================================
Translated from MATLAB model: Models/MP_B/
Reference: Nikravesh - Planar Multibody Dynamics 2nd Edition

Two-body suspension system.  B1 (main link) is constrained with B2 via a
revolute joint at A (J1) and a rev-tran joint (J2, sliding along B1's axis
toward ground O0).  B2 is also pinned to the ground at Q (J3).
A spring connects B1 to O0, and a tire contact force acts at C1 on B1.

Bodies:
  B1: main link  (m=20, J=2.5, r=[0.5840,0.3586], p=6.0819)
  B2: link       (m=2,  J=0.5, r=[0.3450,0.2900], p=0    )

Points (Python 0-based ← MATLAB 1-based):
  [0] A1 : B1 [ 0.00,-0.07]   [1] B1pt: B1 [-0.17, 0.25]
  [2] C1 : B1 [ 0.11,-0.02]   [3] O0  : B0 [ 0.41, 0.83]
  [4] Q0 : B0 [ 0.12, 0.29]   [5] Q2  : B2 [-0.225, 0  ]
  [6] A2 : B2 [ 0.225, 0  ]

Unit Vectors:
  [0] V1: B1, [0.47,-0.88]

Joints:
  j0: rev      A1(B1)   <-> A2(B2)     pt[0] <-> pt[6]
  j1: rev-tran B1pt(B1) ~>  O0(B0)     pt[1] <-> pt[3],  uV[0]
  j2: rev      Q0(B0)   <-> Q2(B2)     pt[4] <-> pt[5]

Forces:
  s0: ptp spring-damper  B1pt <-> O0  (k=20000, L0=0.34, dc=1100)
  s1: user (tire contact at C1 on B1)
  s2: weight
"""

import numpy as np
import os
from PMD.src.builder import *
from PMD.src.solver import *
from PMD.src.mechanics import *
from PMD.examples._plot_utils import plot_comparison

# =============================================================================
# MODEL DEFINITION
# =============================================================================

#%% Bodies
B1 = Body(m=20.0, J=2.5, r=[0.5840, 0.3586], p=6.0819)
B2 = Body(m=2.0,  J=0.5, r=[0.3450, 0.2900], p=0.0)

#%% Points  (MATLAB P1..7 → Python 0..6)
pt_A1  = Point(Bindex=1, sPlocal=np.array([ 0.00, -0.07]))  # [0]
pt_B1  = Point(Bindex=1, sPlocal=np.array([-0.17,  0.25]))  # [1]
pt_C1  = Point(Bindex=1, sPlocal=np.array([ 0.11, -0.02]))  # [2]  ← contact point
pt_O0  = Point(Bindex=0, sPlocal=np.array([ 0.41,  0.83]))  # [3]
pt_Q0  = Point(Bindex=0, sPlocal=np.array([ 0.12,  0.29]))  # [4]
pt_Q2  = Point(Bindex=2, sPlocal=np.array([-0.225, 0.00]))  # [5]
pt_A2  = Point(Bindex=2, sPlocal=np.array([ 0.225, 0.00]))  # [6]

#%% Unit Vectors  (MATLAB V1→[0])
v1 = uVector(Bindex=1, ulocal=np.array([0.47, -0.88]))  # [0]

#%% Joints  (MATLAB J1..3 → Python j0..j2 ; point/uvector indices -1)
j0 = Joint(type='rev',     iPindex=0, jPindex=6)              # A1 <-> A2
j1 = Joint(type='rev-tran', iPindex=1, jPindex=3, iUindex=0)  # B1pt ~> O0 along V1
j2 = Joint(type='rev',     iPindex=4, jPindex=5)              # Q0 <-> Q2

#%% Forces
s0 = Force(type='ptp', iPindex=1, jPindex=3, k=20000.0, L0=0.34, dc=1100.0)

def my_force(B1, s1, pt_C1):
    """Unilateral tire contact — same model as MP_A."""
    del_y = float(B1.r[1]) - s1.L0
    if del_y < 0:
        fy = s1.k * del_y + s1.dc * float(B1.dr[1])
        fsd = np.array([[0.0], [-fy]])
        B1._f += fsd
        B1._n += (s_rot(pt_C1._sP).T @ fsd).item()

s1 = Force(type='user', k=100000.0, L0=0.30, dc=1000.0)
s1.callback = my_force

s2 = Force(type='weight')

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel()
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001),
                    ic_correct=True)

# =============================================================================
# OUTPUT
# =============================================================================

# #%% Save results
# output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', '_test_MP_B.txt')
# nB = model.nB
# nC = model.nC
# nB3 = nB * 3
# header = '\t'.join(['t'] + [f'B{i+1}_{c}' for i in range(nB) for c in ['x', 'y', 'p']])
# np.savetxt(output_file, np.column_stack([T, uT[:, :nB3]]),
#            delimiter='\t', header=header, comments='', fmt='%.8f')
# print(f"[_test_MP_B] Done. nB={nB}, nC={nC}, DOF={nB*3-nC}, points={len(T)}")
# print(f"  Results: {output_file}")

plot_comparison(T, uT, matlab_filename='MP_B.txt', model_title='MP_B')
