"""
MSD (Mass-Spring-Damper) Model - Planar Multibody Dynamics (PMD)
================================================================
KINEMATIC SCHEME & INITIAL STATE
---------------------------------

          +---------+
          |         |
          |    m    |  mass
          |         |
          +---------+
            |     |
            /     |  
            \    _|_
    spring  /   |   | 
            \   |___| damper
            /     |
            |     |
         ===+=====+=== ground

BODY PARAMETERS
---------------
+------+----------+---------------+------------+--------------+
| Body | m  [kg]  |  I  [kg·m²]   | pos  [m]   | vel  [m/s]   |
+------+----------+---------------+------------+--------------+
| mass | 1000.0   |  0.5  (†)     | (0.0, 2.0) | (0.0, 10.0)  |
+------+----------+---------------+------------+--------------+
(†) Placeholder — rotation is blocked by TranJoint, so I has no effect.
"""

import numpy as np
from PMD.src import *
from PMD.gui import PostProcessor


#%% Bodies
mass = Body(
    mass=1000.0,           # kg
    inertia=0.5,           # kg·m² — placeholder, rotation blocked by TranJoint
    position=[0.0, 0.0],   # m  — CM; bottom face at y = 1.5 m = L0
    velocity=[0.0, 10.0],  # m/s — initial velocity upward
    name='mass',
)


#%% Markers
# Ground markers
mk_g_spring = Ground.add_marker([0.0, -2.0])                    # spring lower anchor
mk_g_tran   = Ground.add_marker([0.0, 0.0], theta=np.pi / 2)    # TranJoint guide (Y-axis)

# Mass markers (local frame, CM at origin)
mk_m_spring = mass.add_marker([0.0, -0.5])                    # spring upper anchor (bottom face)
mk_m_tran   = mass.add_marker([0.0,  0.0], theta=np.pi / 2)   # TranJoint guide (Y-axis)


#%% Forces
f_sd = PtpForce(
    iMarker=mk_m_spring,
    jMarker=mk_g_spring,
    k=1.0e6,    # N/m
    L0=1.5,     # m  — natural (undeformed) length
    dc=6500.0,  # N·s/m
    name='spring_damper',
)


#%% Joints
j_gs = TranJoint(
    iMarker=mk_m_tran,
    jMarker=mk_g_tran,
    name='vertical_slider',
)


#%% Simulation
model = PlanarMultibodyModel(
    bodies=[mass],
    joints=[j_gs],
    forces=[f_sd],
)

T, uT = model.solve(
    method='CASADI-DAE',
    t_final=3.0,
    t_eval=np.linspace(0, 3.0, 3001),
    ic_correct=True,
)

post_proc = PostProcessor(model=model, T=T, uT=uT)
post_proc.show()