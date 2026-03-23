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
from PMD.src import *
from PMD.examples._plot_utils import plot_comparison

# =============================================================================
# MODEL DEFINITION
# =============================================================================

#%% Bodies
B1 = Body(mass=20.0, inertia=5.0, position=[0.5, 0.2], orientation=0.0)
B2 = Body(mass=2.0,  inertia=0.5, position=[0.2, 0.1], orientation=0.0)
B3 = Body(mass=2.0,  inertia=0.5, position=[0.8, 0.1], orientation=0.0)

#%% Markers
pt0 = B1.add_marker([-0.3, -0.1])
pt1 = B1.add_marker([ 0.3, -0.1])
pt2 = B2.add_marker([ 0.0,  0.0])
pt3 = B3.add_marker([ 0.0,  0.0])

#%% Functions  (smooth ramp: 0 -> -2*pi rad/s over 2 s)
fn0 = Function(type='c', t_end=2.0, dfdt_end=-2.0*np.pi)

#%% Joints
j0 = RevJoint(iMarker=pt0, jMarker=pt2)
j1 = RevJoint(iMarker=pt1, jMarker=pt3)
j2 = DiscJoint(iBody=B2, R=0.1, x0=0.2)
j3 = DiscJoint(iBody=B3, R=0.1, x0=0.8)
j4 = RelRotJoint(iBody=B2, jBody=B1, iFunct=fn0)

#%% Forces
fw = Weight()

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel(
    bodies=[B1, B2, B3],
    joints=[j0, j1, j2, j3, j4],
    forces=[fw],
    functions=[fn0])
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='Cart_B.txt', model_title='Cart_B')
