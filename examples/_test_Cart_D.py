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

#%% Joints
j0 = Joint(type='rev',  iPoint=pt0, jPoint=pt2)
j1 = Joint(type='rev',  iPoint=pt1, jPoint=pt3)
j2 = Joint(type='disc', iBody=B2, R=0.1, x0=0.2)
j3 = Joint(type='disc', iBody=B3, R=0.1, x0=0.8)

#%% Forces
fw = Force(type='weight')

def my_force():
    """DC motor torque on B2 (reaction on B1) + aerodynamic drag on cart B1."""
    # Motor (DC model capped at T_max)
    omega_max = 4.0 * np.pi
    T_max = 20.0
    omega = abs(B2.dp)
    T_motor = T_max * (1.0 - omega / omega_max)
    if T_motor > T_max:
        T_motor = T_max
    B2._n -= T_motor
    B1._n += T_motor
    # Aerodynamic drag  F = damp_aero * v_x^2  (opposes motion)
    damp_aero = 10.0
    x_d = B1.dr[0, 0]
    f_aero = damp_aero * x_d ** 2
    B1._f[0] -= f_aero   # drag always opposes motion direction

fu = Force(type='user')
fu.callback = my_force

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel(
    bodies=[B1, B2, B3],
    joints=[j0, j1, j2, j3],
    forces=[fw, fu],
    points=[pt0, pt1, pt2, pt3])
T, uT = model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='Cart_D.txt', model_title='Cart_D')
