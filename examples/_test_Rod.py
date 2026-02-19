"""
Rod Model - Planar Multibody Dynamics
======================================
Translated from MATLAB model: Models/Rod/
Reference: Nikravesh - Planar Multibody Dynamics 2nd Edition

Free-flying rigid rod with initial downward velocity, subject to gravity and
damped penalty contact at both endpoints when they reach y=0.

Bodies:
  B1: free rod (m=1, J=0.01, r=[0,1], p=pi/4, dr=[0,-6])

Points:
  pt0: bottom endpoint, B1-local [0,-1]
  pt1: top endpoint,    B1-local [0, 1]

Joints:
  (none - free body)

Forces:
  f0: weight (gravity)
  f1: user (damped penalty contact at rod endpoints)

Contact model: F_n = k_c * depth + c_c * depth_rate
  k_c = 1e4   (contact stiffness)
  c_c = 200   (contact damping; ~ ζ=1 critical damping for k_c, m=1)
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
B1 = Body(m=1.0, J=0.01, r=[0.0, 1.0], p=np.pi/4, dr=[0.0, -6.0])

#%% Forces
f0 = Force(type='weight')   # gravity

_k_c = 1e4   # contact stiffness
_c_c = 200.0 # contact damping (critically damped: 2*sqrt(k_c*m))

def my_force(B1):
    """Damped penalty contact at both rod endpoints (ground at y=0).

    Velocity of endpoint P on a rigid body:
        v_P_y = v_CM_y + omega * sp_glob_x      (2D rigid body formula)
    where sp_glob = A(phi) @ sp_local.
    """
    A_phi  = A_matrix(float(B1.p))
    omega  = float(B1.dp)
    vy_cm  = float(B1.dr[1])
    for sp_local in [np.array([0., -1.]), np.array([0., 1.])]:
        sp_glob = A_phi @ sp_local                    # global offset from CM
        p_y      = float(B1.r[1]) + sp_glob[1]        # global y of endpoint
        if p_y < 0:
            v_y  = vy_cm + omega * sp_glob[0]         # vertical vel of endpoint
            depth = -p_y
            Fn    = _k_c * depth - _c_c * v_y         # upward contact force
            if Fn < 0:
                Fn = 0.0                               # no tensile contact
            B1._f[1] += Fn
            B1._n    += float(sp_glob[0]) * Fn         # torque: r_x * F_y

f1 = Force(type='user')
f1.callback = my_force

# =============================================================================
# SIMULATION
# =============================================================================

#%% Create model and solve
model = PlanarMultibodyModel()
T, uT = model.solve(method='Radau', t_final=3.0, t_eval=np.linspace(0, 3, 601))

# =============================================================================
# OUTPUT
# =============================================================================

#%% Save results
output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_test_Rod.txt')
nB = model.nB
nC = model.nC
nB3 = nB * 3
header = '\t'.join(['t'] + [f'B{i+1}_{c}' for i in range(nB) for c in ['x', 'y', 'p']])
np.savetxt(output_file, np.column_stack([T, uT[:, :nB3]]),
           delimiter='\t', header=header, comments='', fmt='%.8f')
print(f"[_test_Rod] Done. nB={nB}, nC={nC}, DOF={nB*3-nC}, points={len(T)}")
print(f"  Results: {output_file}")
