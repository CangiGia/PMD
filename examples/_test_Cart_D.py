"""
Cart_D Model - Planar Multibody Dynamics
=========================================
Translated from MATLAB model: Models/Cart_D/
Reference: Nikravesh - Planar Multibody Dynamics 2nd Edition

Same as Cart_C (motor driven) but with an added aerodynamic drag force
on the cart body B1 proportional to the square of its velocity.

Bodies, Points, Joints: identical to Cart_C.

Forces:
  fw: weight
  fu: user (DC motor torque + aerodynamic drag)
"""

import numpy as np
import os
from pmd.core import *
from examples._plot_utils import plot_comparison

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

#%% Joints
j0 = RevJoint(iMarker=pt0, jMarker=pt2)
j1 = RevJoint(iMarker=pt1, jMarker=pt3)
j2 = DiscJoint(iBody=B2, R=0.1, x0=0.2)
j3 = DiscJoint(iBody=B3, R=0.1, x0=0.8)

#%% Forces
fw = Weight()

def dc_motor_and_drag():
    """DC motor torque on B2 (reaction on B1) + aerodynamic drag on cart B1."""
    # Motor (DC model capped at T_max)
    omega_max = 4.0 * np.pi
    T_max = 20.0
    omega = abs(B2.angular_velocity)
    T_motor = T_max * (1.0 - omega / omega_max)
    if T_motor > T_max:
        T_motor = T_max
    # Aerodynamic drag  F = damp_aero * v_x^2  (opposes motion)
    damp_aero = 10.0
    x_d = B1.velocity[0, 0]
    f_aero = damp_aero * x_d ** 2
    return [
        {'body': B2, 'force': [0, 0], 'torque': -T_motor},
        {'body': B1, 'force': [-f_aero, 0], 'torque': T_motor},
    ]

fu = UserForce(callback=dc_motor_and_drag)

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel(
    bodies=[B1, B2, B3],
    joints=[j0, j1, j2, j3],
    forces=[fw, fu])
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='Cart_D.txt', model_title='Cart_D')
