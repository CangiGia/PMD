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

#%% Joints
j0 = Joint(type='rev',  iPindex=0, jPindex=2)
j1 = Joint(type='rev',  iPindex=1, jPindex=3)
j2 = Joint(type='disc', iBindex=2, R=0.1, x0=0.2)
j3 = Joint(type='disc', iBindex=3, R=0.1, x0=0.8)

#%% Forces
fw = Force(type='weight')

def my_force(B1, B2):
    """DC motor torque on B2 (reaction on B1) + aerodynamic drag on cart B1."""
    # Motor (DC model capped at T_max)
    omega_max = 4.0 * np.pi
    T_max = 20.0
    omega = abs(float(B2.dp))
    T_motor = T_max * (1.0 - omega / omega_max)
    if T_motor > T_max:
        T_motor = T_max
    B2._n -= T_motor
    B1._n += T_motor
    # Aerodynamic drag  F = damp_aero * v_x^2  (opposes motion)
    damp_aero = 10.0
    x_d = float(B1.dr[0])
    f_aero = damp_aero * x_d ** 2
    B1._f[0] -= f_aero   # drag always opposes motion direction

fu = Force(type='user')
fu.callback = my_force

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
output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_test_Cart_D.txt')
nB = model.nB
nC = model.nC
nB3 = nB * 3
header = '\t'.join(['t'] + [f'B{i+1}_{c}' for i in range(nB) for c in ['x', 'y', 'p']])
np.savetxt(output_file, np.column_stack([T, uT[:, :nB3]]),
           delimiter='\t', header=header, comments='', fmt='%.8f')
print(f"[_test_Cart_D] Done. nB={nB}, nC={nC}, DOF={nB*3-nC}, points={len(T)}")
print(f"  Results: {output_file}")
