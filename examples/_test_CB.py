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

# =============================================================================
# MODEL DEFINITION
# =============================================================================

#%% Bodies
B1 = Body(m=1.0, J=1.0, r=[1.0, 0.2], p=0.0, dr=[0.0, 0.0])

#%% Points  (MATLAB P1→pt0, P2→pt1  ;  MATLAB iPindex=2→Python ipindex=1, jPindex=1→0)
pt0 = Point(Bindex=0, sPlocal=np.array([0.0, 0.2]))  # wall anchor on ground
pt1 = Point(Bindex=1, sPlocal=np.array([0.0, 0.0]))  # body reference on B1

#%% Unit Vectors  (MATLAB U1→u0, U2→u1  ; MATLAB iUindex=2→Python 1, jUindex=1→0)
u0 = uVector(Bindex=0, ulocal=np.array([1.0, 0.0]))  # horizontal ground vector
u1 = uVector(Bindex=1, ulocal=np.array([1.0, 0.0]))  # horizontal body vector

#%% Joints
j0 = Joint(type='tran', iPindex=1, jPindex=0, iUindex=1, jUindex=0)

#%% Forces
f0 = Force(type='ptp', iPindex=1, jPindex=0, k=10.0, L0=0.8, dc=0.0)
f1 = Force(type='weight')

def my_force(B1):
    """Anderson friction model — body slides against conveyor belt at v_conv=0.1 m/s."""
    mu_d = 0.15; mu_s = 0.2; mu_v = 0.0
    v_s = 0.001; p = 2; k_t = 10000
    fy_normal = 9.81          # normal force (body weight on horizontal surface)
    v_conv = 0.1
    v = v_conv - float(B1.dr[0])   # relative sliding velocity
    ff = friction_A(mu_s, mu_d, v_s, p, k_t, v, fy_normal)
    fx = ff + mu_v * v * fy_normal
    B1._f[0] += fx

f2 = Force(type='user')
f2.callback = my_force

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel()
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

# =============================================================================
# OUTPUT
# =============================================================================

#%% Save results
output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', '_test_CB.txt')
nB = model.nB
nC = model.nC
nB3 = nB * 3
header = '\t'.join(['t'] + [f'B{i+1}_{c}' for i in range(nB) for c in ['x', 'y', 'p']])
np.savetxt(output_file, np.column_stack([T, uT[:, :nB3]]),
           delimiter='\t', header=header, comments='', fmt='%.8f')
print(f"[_test_CB] Done. nB={nB}, nC={nC}, DOF={nB*3-nC}, points={len(T)}")
print(f"  Results: {output_file}")
