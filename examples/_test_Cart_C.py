"""
Cart_C Model - Planar Multibody Dynamics
=========================================
Translated from MATLAB model: Models/Cart_C/
Reference: Nikravesh - Planar Multibody Dynamics 2nd Edition

Cart on disc wheels driven by a user-defined electric motor (DC motor model).
The motor applies torque to wheel B2 with reaction on cart B1.  No rel-rot
joint is used; the motor is purely a user force.

Bodies:
  B1: cart  (m=20, J=5,   r=[0.5, 0.2], p=0)
  B2: left  (m=2,  J=0.5, r=[0.2, 0.1], p=0)
  B3: right (m=2,  J=0.5, r=[0.8, 0.1], p=0)

Joints:
  j0: rev   pt0 <-> pt2
  j1: rev   pt1 <-> pt3
  j2: disc  B2, R=0.1, x0=0.2
  j3: disc  B3, R=0.1, x0=0.8

Forces:
  fw: weight
  fu: user (DC motor torque on B2 / reaction on B1)
"""

import numpy as np
import os
from PMD.src import *
from PMD.examples._plot_utils import plot_comparison

# =============================================================================
# MODEL DEFINITION
# =============================================================================

#%% Bodies
B1 = Body(m=20.0, J=5.0, r=[0.5, 0.2], p=0.0)
B2 = Body(m=2.0,  J=0.5, r=[0.2, 0.1], p=0.0)
B3 = Body(m=2.0,  J=0.5, r=[0.8, 0.1], p=0.0)

#%% Markers
pt0 = B1.add_marker([-0.3, -0.1])
pt1 = B1.add_marker([ 0.3, -0.1])
pt2 = B2.add_marker([ 0.0,  0.0])
pt3 = B3.add_marker([ 0.0,  0.0])

#%% Joints (4 joints -- no rel-rot)
j0 = Joint(type='rev',  iMarker=pt0, jMarker=pt2)
j1 = Joint(type='rev',  iMarker=pt1, jMarker=pt3)
j2 = Joint(type='disc', iBody=B2, R=0.1, x0=0.2)
j3 = Joint(type='disc', iBody=B3, R=0.1, x0=0.8)

#%% Forces
fw = Force(type='weight')

def my_force():
    """DC motor: torque on B2 with equal and opposite reaction on B1."""
    omega_max = 4.0 * np.pi   # [rad/s] no-load speed
    T_max = 20.0               # [N*m]  stall torque
    omega = abs(float(B2.dp))
    T_motor = T_max * (1.0 - omega / omega_max)
    if T_motor > T_max:
        T_motor = T_max
    B2._n -= T_motor
    B1._n += T_motor

fu = Force(type='user')
fu.callback = my_force

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
    plot_comparison(T, uT, matlab_filename='Cart_C.txt', model_title='Cart_C')
