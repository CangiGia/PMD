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
from pmd.core import *
from examples._plot_utils import plot_comparison

# =============================================================================
# MODEL DEFINITION
# =============================================================================

#%% Bodies
B1 = Body(mass=20.0, inertia=2.5, position=[0.5840, 0.3586], orientation=6.0819)

#%% Markers
pt_A1  = B1.add_marker([ 0.00, -0.07])                          # [0]
pt_B1  = B1.add_marker([-0.17,  0.25], theta=np.arctan2(-0.88, 0.47))  # [1] + V1 orientation
pt_C1  = B1.add_marker([ 0.11, -0.02])                          # [2] <- contact point
pt_O0  = Ground.add_marker([ 0.41,  0.83])                      # [3]
pt_Q0  = Ground.add_marker([ 0.12,  0.29])                      # [4]

#%% Joints  (MATLAB J1..2 -> Python j0..j1)
j0 = RevRevJoint(iMarker=pt_A1, jMarker=pt_Q0, L=0.45)     # A1 <~~> Q0
j1 = RevTranJoint(iMarker=pt_B1, jMarker=pt_O0)              # B1pt ~> O0 along V1

#%% Forces
s0 = PtpForce(iMarker=pt_B1, jMarker=pt_O0, k=20000.0, L0=0.34, dc=1100.0)

_k_tire = 100000.0
_L0_tire = 0.30
_dc_tire = 1000.0

def tire_contact():
    """Unilateral tire contact -- same model as MP_A/B."""
    del_y = B1.position[1, 0] - _L0_tire
    if del_y < 0:
        fy = _k_tire * del_y + _dc_tire * B1.velocity[1, 0]
        fsd = np.array([[0.0], [-fy]])
        torque = (rotate_90(pt_C1._sP).T @ fsd).item()
        return [{'body': B1, 'force': [0.0, -fy], 'torque': torque}]
    return []

s1 = UserForce(callback=tire_contact)

s2 = Weight()

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel(
    bodies=[B1],
    joints=[j0, j1],
    forces=[s0, s1, s2])
T, uT = model.solve(t_final=10.0, t_eval=np.linspace(0, 10, 10001),
                    ic_correct=True)

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='MP_C.txt', model_title='MP_C')
