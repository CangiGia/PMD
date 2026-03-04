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
from PMD.examples._plot_utils import plot_comparison

# =============================================================================
# MODEL DEFINITION
# =============================================================================

#%% Bodies
B1 = Body(m=20.0, J=5.0, r=[0.5, 0.2], p=0.0)
B2 = Body(m=2.0,  J=0.5, r=[0.2, 0.1], p=0.0)
B3 = Body(m=2.0,  J=0.5, r=[0.8, 0.1], p=0.0)

#%% Points
pt0 = Point(body=B1, sPlocal=np.array([-0.3, -0.1]))
pt1 = Point(body=B1, sPlocal=np.array([ 0.3, -0.1]))
pt2 = Point(body=B2, sPlocal=np.array([ 0.0,  0.0]))
pt3 = Point(body=B3, sPlocal=np.array([ 0.0,  0.0]))

#%% Functions  (smooth ramp: 0 -> -2*pi rad/s over 2 s)
fn0 = Function(type='c', t_end=2.0, dfdt_end=-2.0*np.pi)

#%% Joints
j0 = Joint(type='rev',     iPoint=pt0, jPoint=pt2)
j1 = Joint(type='rev',     iPoint=pt1, jPoint=pt3)
j2 = Joint(type='disc',    iBody=B2, R=0.1, x0=0.2)
j3 = Joint(type='disc',    iBody=B3, R=0.1, x0=0.8)
j4 = Joint(type='rel-rot', iBody=B2, jBody=B1, iFunct=fn0)

#%% Forces
fw = Force(type='weight')

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel(
    bodies=[B1, B2, B3],
    joints=[j0, j1, j2, j3, j4],
    forces=[fw],
    points=[pt0, pt1, pt2, pt3],
    functions=[fn0])
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='Cart_B.txt', model_title='Cart_B')
