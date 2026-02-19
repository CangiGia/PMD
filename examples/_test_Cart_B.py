"""
Cart_B Model - Planar Multibody Dynamics
=========================================
Translated from MATLAB model: Models/Cart_B/
Reference: Nikravesh - Planar Multibody Dynamics 2nd Edition

Same cart-on-wheels as Cart_A, but the motor uses a smooth cosine-based
ramp-up (type-'c' function): it starts from rest and reaches -2*pi rad/s
over 2 seconds, then runs at constant speed.

Bodies, Points, and Joints: identical to Cart_A.

Functions:
  fn0: type='c', t_end=2.0, dfdt_end=-2*pi  (smooth velocity ramp-up)

Forces:
  fw: weight
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
B1 = Body(m=20.0, J=5.0, r=[0.5, 0.2], p=0.0)
B2 = Body(m=2.0,  J=0.5, r=[0.2, 0.1], p=0.0)
B3 = Body(m=2.0,  J=0.5, r=[0.8, 0.1], p=0.0)

#%% Points
pt0 = Point(Bindex=1, sPlocal=np.array([-0.3, -0.1]))
pt1 = Point(Bindex=1, sPlocal=np.array([ 0.3, -0.1]))
pt2 = Point(Bindex=2, sPlocal=np.array([ 0.0,  0.0]))
pt3 = Point(Bindex=3, sPlocal=np.array([ 0.0,  0.0]))

#%% Functions  (smooth ramp: 0 → -2*pi rad/s over 2 s)
fn0 = Function(type='c', t_end=2.0, dfdt_end=-2.0*np.pi)

#%% Joints
j0 = Joint(type='rev',     iPindex=0, jPindex=2)
j1 = Joint(type='rev',     iPindex=1, jPindex=3)
j2 = Joint(type='disc',    iBindex=2, R=0.1, x0=0.2)
j3 = Joint(type='disc',    iBindex=3, R=0.1, x0=0.8)
j4 = Joint(type='rel-rot', iBindex=2, jBindex=1, iFunct=1)

#%% Forces
fw = Force(type='weight')

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel()
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 1001))

# =============================================================================
# OUTPUT
# =============================================================================

#%% Save results
output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_test_Cart_B.txt')
nB = model.nB
nC = model.nC
nB3 = nB * 3
header = '\t'.join(['t'] + [f'B{i+1}_{c}' for i in range(nB) for c in ['x', 'y', 'p']])
np.savetxt(output_file, np.column_stack([T, uT[:, :nB3]]),
           delimiter='\t', header=header, comments='', fmt='%.8f')
print(f"[_test_Cart_B] Done. nB={nB}, nC={nC}, DOF={nB*3-nC}, points={len(T)}")
print(f"  Results: {output_file}")
