"""
MP_C Model - Planar Multibody Dynamics
========================================
Translated from MATLAB model: Models/MP_C/
Reference: Nikravesh - Planar Multibody Dynamics 2nd Edition

Single-body suspension system.  Body B1 is constrained to the ground via:
  - a rev-rev joint (constant-length link Q0-A1, J1)
  - a rev-tran joint (B1pt slides along its own axis toward O0, J2)
A spring connects B1 to O0, and a tire contact force acts at C1.

Bodies:
  B1: main link  (m=20, J=2.5, r=[0.5840,0.3586], p=6.0819)

Points (Python 0-based <- MATLAB 1-based):
  [0] A1 : B1 [ 0.00,-0.07]   [1] B1pt: B1 [-0.17, 0.25]
  [2] C1 : B1 [ 0.11,-0.02]   [3] O0  : B0 [ 0.41, 0.83]
  [4] Q0 : B0 [ 0.12, 0.29]

Unit Vectors:
  [0] V1: B1, [0.47,-0.88]

Joints:
  j0: rev-rev   A1(B1) <~~ Q0(B0)    pt[0] <-> pt[4],  L=0.45
  j1: rev-tran  B1pt(B1)~> O0(B0)    pt[1] <-> pt[3],  uV[0]

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

#%% Points  (MATLAB P1..5 -> Python 0..4)
pt_A1  = Point(body=B1, sPlocal=np.array([ 0.00, -0.07]))  # [0]
pt_B1  = Point(body=B1, sPlocal=np.array([-0.17,  0.25]))  # [1]
pt_C1  = Point(body=B1, sPlocal=np.array([ 0.11, -0.02]))  # [2]  <- contact point
pt_O0  = Point(body=Ground, sPlocal=np.array([ 0.41,  0.83]))  # [3]
pt_Q0  = Point(body=Ground, sPlocal=np.array([ 0.12,  0.29]))  # [4]

#%% Unit Vectors  (MATLAB V1->[0])
v1 = uVector(body=B1, ulocal=np.array([0.47, -0.88]))  # [0]

#%% Joints  (MATLAB J1..2 -> Python j0..j1)
j0 = Joint(type='rev-rev',  iPoint=pt_A1, jPoint=pt_Q0, L=0.45)     # A1 <~~> Q0
j1 = Joint(type='rev-tran', iPoint=pt_B1, jPoint=pt_O0, iUvec=v1)   # B1pt ~> O0 along V1

#%% Forces
s0 = Force(type='ptp', iPoint=pt_B1, jPoint=pt_O0, k=20000.0, L0=0.34, dc=1100.0)

s1 = Force(type='user', k=100000.0, L0=0.30, dc=1000.0)

def my_force():
    """Unilateral tire contact -- same model as MP_A/B."""
    del_y = B1.r[1, 0] - s1.L0
    if del_y < 0:
        fy = s1.k * del_y + s1.dc * B1.dr[1, 0]
        fsd = np.array([[0.0], [-fy]])
        B1._f += fsd
        B1._n += (s_rot(pt_C1._sP).T @ fsd).item()

s1.callback = my_force

s2 = Force(type='weight')

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel(
    bodies=[B1],
    joints=[j0, j1],
    forces=[s0, s1, s2],
    points=[pt_A1, pt_B1, pt_C1, pt_O0, pt_Q0],
    uvectors=[v1])
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001),
                    ic_correct=True)

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='MP_C.txt', model_title='MP_C')
