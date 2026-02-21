"""
Cart_A Model - Planar Multibody Dynamics
=========================================
Translated from MATLAB model: Models/Cart_A/
Reference: Nikravesh - Planar Multibody Dynamics 2nd Edition

Cart supported on two disc wheels. A motor drives wheel B2 relative to the
cart B1 at a constant angular velocity (-2*pi rad/s, type-'a' polynomial).
The wheel-ground contact is enforced via disc constraints.


Bodies:
  B1: cart body  (m=20, J=5,   r=[0.5, 0.2],  p=0)
  B2: left wheel (m=2,  J=0.5, r=[0.2, 0.1],  p=0)
  B3: right wheel(m=2,  J=0.5, r=[0.8, 0.1],  p=0)

Points:
  pt0: B1 left  axle,  B1-local [-0.3,-0.1]
  pt1: B1 right axle,  B1-local [ 0.3,-0.1]
  pt2: B2 center,      B2-local [ 0,   0  ]
  pt3: B3 center,      B3-local [ 0,   0  ]

Joints:
  j0: rev   B1-left  <-> B2-center  (pt0 <-> pt2)
  j1: rev   B1-right <-> B3-center  (pt1 <-> pt3)
  j2: disc  B2, R=0.1, x0=0.2
  j3: disc  B3, R=0.1, x0=0.8
  j4: rel-rot  B2 relative to B1, driven by fn0

Functions:
  fn0: type='a', f(t)=-2*pi*t  (constant angular velocity motor)

Forces:
  fw: weight
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
B1 = Body(m=20.0, J=5.0,  r=[0.5, 0.2], p=0.0)
B2 = Body(m=2.0,  J=0.5,  r=[0.2, 0.1], p=0.0)
B3 = Body(m=2.0,  J=0.5,  r=[0.8, 0.1], p=0.0)

#%% Points  (MATLAB P1..4 → Python 0..3)
pt0 = Point(Bindex=1, sPlocal=np.array([-0.3, -0.1]))  # B1 left  axle
pt1 = Point(Bindex=1, sPlocal=np.array([ 0.3, -0.1]))  # B1 right axle
pt2 = Point(Bindex=2, sPlocal=np.array([ 0.0,  0.0]))  # B2 center
pt3 = Point(Bindex=3, sPlocal=np.array([ 0.0,  0.0]))  # B3 center

#%% Functions  (MATLAB Functs(1) → Python fn0, iFunct=1)
fn0 = Function(type='a', coeff=[0.0, -2.0*np.pi, 0.0])

#%% Joints  (MATLAB J1..5 → Python j0..j4)
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
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001),
                    ic_correct=True)

# =============================================================================
# OUTPUT
# =============================================================================

# #%% Save results
# output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', '_test_Cart_A.txt')
# nB = model.nB
# nC = model.nC
# nB3 = nB * 3
# header = '\t'.join(['t'] + [f'B{i+1}_{c}' for i in range(nB) for c in ['x', 'y', 'p']])
# np.savetxt(output_file, np.column_stack([T, uT[:, :nB3]]),
#            delimiter='\t', header=header, comments='', fmt='%.8f')
# print(f"[_test_Cart_A] Done. nB={nB}, nC={nC}, DOF={nB*3-nC}, points={len(T)}")
# print(f"  Results: {output_file}")

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='Cart_A.txt', model_title='Cart_A')
