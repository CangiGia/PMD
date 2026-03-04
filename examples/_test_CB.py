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
from PMD.src.builder import *
from PMD.src.solver import *
from PMD.src.mechanics import *
from PMD.examples._plot_utils import plot_comparison

# =============================================================================
# MODEL DEFINITION
# =============================================================================

#%% Bodies
B1 = Body(m=1.0, J=1.0, r=[1.0, 0.2], p=0.0, dr=[0.0, 0.0])

#%% Points
pt0 = Point(body=Ground, sPlocal=np.array([0.0, 0.2]))  # wall anchor on ground
pt1 = Point(body=B1, sPlocal=np.array([0.0, 0.0]))  # body reference on B1

#%% Unit Vectors
u0 = uVector(body=Ground, ulocal=np.array([1.0, 0.0]))  # horizontal ground vector
u1 = uVector(body=B1, ulocal=np.array([1.0, 0.0]))  # horizontal body vector

#%% Joints
j0 = Joint(type='tran', iPoint=pt1, jPoint=pt0, iUvec=u1, jUvec=u0)

#%% Forces
f0 = Force(type='ptp', iPoint=pt1, jPoint=pt0, k=10.0, L0=0.8, dc=0.0)
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
    forces=[f0, f1, f2],
    points=[pt0, pt1],
    uvectors=[u0, u1])
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='CB.txt', model_title='CB')
