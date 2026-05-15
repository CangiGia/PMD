r"""
MSD (Mass-Spring-Damper) Model - Planar Multibody Dynamics (pmd)
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
    spring  /   |___| 
            \   | | |damper
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
from pmd.core import (
    Body,
    Ground,
    PlanarMultibodyModel,
    PtpForce,
    TranJoint,
)
from pmd.core.shapes import Rectangle
from pmd.core.model import _GroundType


T_FINAL = 3.0
N_EVAL = 3001
IC_CORRECT = True

if _GroundType._instance is not None:
    _GroundType._markers = [_GroundType._instance.origin]

mass = Body(
    mass=1000.0,
    inertia=0.5,           # placeholder — rotation blocked by TranJoint
    position=[0.0, 0.0],
    velocity=[0.0, 10.0],
    name='mass',
    shape=Rectangle(width=0.4, height=0.3),
)

mk_g_spring = Ground.add_marker([0.0, -2.0])
mk_g_tran   = Ground.add_marker([0.0, 0.0], theta=np.pi / 2)
mk_m_spring = mass.add_marker([0.0, -0.5])
mk_m_tran   = mass.add_marker([0.0,  0.0], theta=np.pi / 2)

f_sd = PtpForce(
    iMarker=mk_m_spring,
    jMarker=mk_g_spring,
    k=1.0e6,
    L0=1.5,     # natural (undeformed) length [m]
    dc=6500.0,
    name='spring_damper',
)

j_gs = TranJoint(
    iMarker=mk_m_tran,
    jMarker=mk_g_tran,
    name='vertical_slider',
)

model = PlanarMultibodyModel(
    bodies=[mass],
    joints=[j_gs],
    forces=[f_sd],
)


if __name__ == "__main__":
    from pmd.gui import PostProcessor

    T, uT = model.solve(
        t_final=T_FINAL,
        t_eval=np.linspace(0, T_FINAL, N_EVAL),
        ic_correct=IC_CORRECT,
    )
    PostProcessor(model=model, T=T, uT=uT).show()
