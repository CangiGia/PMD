"""
MP_A Model - Planar Multibody Dynamics
========================================
Translated from MATLAB model: Models/MP_A/
Reference: Nikravesh - Planar Multibody Dynamics 2nd Edition

Three-body planar suspension system.  Body B1 (main suspension link) is
connected to the ground via a revolute joint at A (J1) and to a link B2
via revolute joints at Q (J2) and O (J3).  A translational joint (J4)
constrains B1 relative to B3 (second link).  A spring connects B1 to the
fixed ground point O0.  A tire contact force acts at point C1 on B1.

Initial conditions are prescribed at the MATLAB-validated equilibrium.

Bodies:
  B1: main body  (m=20, J=2.5, r=[0.5840,0.3586], p=6.0819)
  B2: link       (m=2,  J=0.5, r=[0.3450,0.2900], p=0    )
  B3: pendulum   (m=0.5,J=0.2, r=[0.4528,0.6862], p=5.0019)

Points (Python 0-based <- MATLAB 1-based):
  [0] A1 : B1 [-0.00,-0.07]   [1] B1pt: B1 [-0.17, 0.25]
  [2] C1 : B1 [ 0.11,-0.02]   [3] O0  : B0 [ 0.41, 0.83]
  [4] Q0 : B0 [ 0.12, 0.29]   [5] Q2  : B2 [-0.225, 0  ]
  [6] A2 : B2 [ 0.225, 0  ]   [7] O3  : B3 [-0.15,  0  ]

Unit Vectors:
  [0] V1: B1, [0.47,-0.88]
  [1] V2: B3, [1, 0]

Joints:
  j0: rev   A1(B1) <-> A2(B2)     pt[0]  <-> pt[6]
  j1: rev   Q0(B0) <-> Q2(B2)     pt[4]  <-> pt[5]
  j2: rev   O0(B0) <-> O3(B3)     pt[3]  <-> pt[7]
  j3: tran  B1pt(B1)<->O3(B3)     pt[1]  <-> pt[7],  uV[0] // uV[1]

Forces:
  s0: ptp spring-damper  B1pt<->O0  (k=20000, L0=0.34, dc=1100)
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
B3 = Body(m=0.5,  J=0.2, r=[0.4528, 0.6862], p=5.0019)

#%% Points  (MATLAB P1..8 -> Python 0..7)
pt_A1  = Point(body=B1, sPlocal=np.array([ 0.00, -0.07]))  # [0]
pt_B1  = Point(body=B1, sPlocal=np.array([-0.17,  0.25]))  # [1]
pt_C1  = Point(body=B1, sPlocal=np.array([ 0.11, -0.02]))  # [2]  <- contact point
pt_O0  = Point(body=Ground, sPlocal=np.array([ 0.41,  0.83]))  # [3]
pt_Q0  = Point(body=Ground, sPlocal=np.array([ 0.12,  0.29]))  # [4]
pt_Q2  = Point(body=B2, sPlocal=np.array([-0.225, 0.00]))  # [5]
pt_A2  = Point(body=B2, sPlocal=np.array([ 0.225, 0.00]))  # [6]
pt_O3  = Point(body=B3, sPlocal=np.array([-0.15,  0.00]))  # [7]

#%% Unit Vectors  (MATLAB V1->[0], V2->[1])
v1 = uVector(body=B1, ulocal=np.array([0.47, -0.88]))  # [0]  body 1 direction
v2 = uVector(body=B3, ulocal=np.array([1.0,   0.0]))   # [1]  body 3 horizontal

#%% Joints  (MATLAB J1..4 -> Python j0..j3)
j0 = Joint(type='rev',  iPoint=pt_A1, jPoint=pt_A2)                         # A1 <-> A2
j1 = Joint(type='rev',  iPoint=pt_Q0, jPoint=pt_Q2)                         # Q0 <-> Q2
j2 = Joint(type='rev',  iPoint=pt_O0, jPoint=pt_O3)                         # O0 <-> O3
j3 = Joint(type='tran', iPoint=pt_B1, jPoint=pt_O3, iUvec=v1, jUvec=v2)     # B1pt tran O3

#%% Forces
s0 = Force(type='ptp', iPoint=pt_B1, jPoint=pt_O0, k=20000.0, L0=0.34, dc=1100.0)

s1 = Force(type='user', k=100000.0, L0=0.30, dc=1000.0)

def my_force():
    """Unilateral tire contact: penalty spring+damper active when B1.y < contact radius."""
    del_y = B1.r[1, 0] - s1.L0   # y_B1 - contact_radius
    if del_y < 0:
        fy = s1.k * del_y + s1.dc * B1.dr[1, 0]
        fsd = np.array([[0.0], [-fy]])      # contact force (upward)
        B1._f += fsd
        B1._n += (s_rot(pt_C1._sP).T @ fsd).item()

s1.callback = my_force

s2 = Force(type='weight')

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel(
    bodies=[B1, B2, B3],
    joints=[j0, j1, j2, j3],
    forces=[s0, s1, s2],
    points=[pt_A1, pt_B1, pt_C1, pt_O0, pt_Q0, pt_Q2, pt_A2, pt_O3],
    uvectors=[v1, v2])
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001),
                    ic_correct=True)

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='MP_A.txt', model_title='MP_A')
