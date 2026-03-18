"""
CB (Sliding Body) Model - Planar Multibody Dynamics
====================================================
Translated from MATLAB model: Models/CB/
Reference: Nikravesh - Planar Multibody Dynamics 2nd Edition

A single rigid body sliding along a horizontal rail, connected to a fixed
wall by a spring, subject to Anderson friction model.  The body tries to
maintain a conveyor-belt velocity v_conv = 0.1 m/s.

Bodies:
  B1: slider (m=1, J=1, r=[1, 0.2], p=0)

Points:
  pt0: wall anchor,  ground-fixed  [0, 0.2]
  pt1: body center,  B1-fixed      [0, 0]

Unit Vectors:
  u0: ground horizontal  [1, 0]
  u1: body horizontal    [1, 0]

Joints:
  j0: translational  (pt1 <-> pt0,  u1 <-> u0)

Forces:
  f0: ptp spring-damper  (pt1-pt0, k=10, L0=0.8, dc=0)
  f1: weight
  f2: user  (Anderson friction model)
"""

import numpy as np
import os
from PMD.src import *
from PMD.examples._plot_utils import plot_comparison

# =============================================================================
# MODEL DEFINITION
# =============================================================================

#%% Bodies
B1 = Body(m=1.0, J=1.0, r=[1.0, 0.2], p=0.0, dr=[0.0, 0.0])

#%% Markers (merge Point + uVector into single Marker with theta for tran joint)
# pt0 was Point(body=Ground, [0,0.2]) and u0 was uVector(body=Ground, [1,0]) -> theta=0.0
pt0 = Ground.add_marker([0.0, 0.2], theta=0.0)  # wall anchor + horizontal direction on ground

# pt1 was Point(body=B1, [0,0]) and u1 was uVector(body=B1, [1,0]) -> theta=0.0
pt1 = B1.add_marker([0.0, 0.0], theta=0.0)  # body reference + horizontal direction on B1

#%% Joints
j0 = Joint(type='tran', iMarker=pt1, jMarker=pt0)

#%% Forces
f0 = Force(type='ptp', iMarker=pt1, jMarker=pt0, k=10.0, L0=0.8, dc=0.0)
f1 = Force(type='weight')

def my_force():
    """Anderson friction model -- body slides against conveyor belt at v_conv=0.1 m/s."""
    mu_d = 0.15; mu_s = 0.2; mu_v = 0.0
    v_s = 0.001; p = 2; k_t = 10000
    fy_normal = 9.81          # normal force (body weight on horizontal surface)
    v_conv = 0.1
    v = v_conv - B1.dr[0, 0]   # relative sliding velocity
    ff = friction_A(mu_s, mu_d, v_s, p, k_t, v, fy_normal)
    fx = ff + mu_v * v * fy_normal
    B1._f[0] += fx

f2 = Force(type='user')
f2.callback = my_force

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel(
    bodies=[B1],
    joints=[j0],
    forces=[f0, f1, f2])
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='CB.txt', model_title='CB')
