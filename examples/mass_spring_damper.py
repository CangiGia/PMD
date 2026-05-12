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
from pmd.core import (
    Body,
    Ground,
    PlanarMultibodyModel,
    PtpForce,
    TranJoint,
)
from pmd.core.model import _GroundType


T_FINAL = 3.0
N_EVAL = 3001
IC_CORRECT = True


def build_model() -> PlanarMultibodyModel:
    """Build the mass-spring-damper :class:`PlanarMultibodyModel` (no solve)."""
    # Reset Ground singleton so repeated builds in the same process are independent.
    gt = _GroundType._instance
    if gt is not None:
        _GroundType._markers = [gt.origin]

    # Bodies
    mass = Body(
        mass=1000.0,           # kg
        inertia=0.5,           # kg·m² — placeholder, rotation blocked by TranJoint
        position=[0.0, 0.0],   # m  — CM; bottom face at y = 1.5 m = L0
        velocity=[0.0, 10.0],  # m/s — initial velocity upward
        name='mass',
    )

    # Markers
    mk_g_spring = Ground.add_marker([0.0, -2.0])                    # spring lower anchor
    mk_g_tran   = Ground.add_marker([0.0, 0.0], theta=np.pi / 2)    # TranJoint guide (Y-axis)
    mk_m_spring = mass.add_marker([0.0, -0.5])                      # spring upper anchor
    mk_m_tran   = mass.add_marker([0.0,  0.0], theta=np.pi / 2)     # TranJoint guide (Y-axis)

    # Forces
    f_sd = PtpForce(
        iMarker=mk_m_spring,
        jMarker=mk_g_spring,
        k=1.0e6,    # N/m
        L0=1.5,     # m  — natural (undeformed) length
        dc=6500.0,  # N·s/m
        name='spring_damper',
    )

    # Joints
    j_gs = TranJoint(
        iMarker=mk_m_tran,
        jMarker=mk_g_tran,
        name='vertical_slider',
    )

    return PlanarMultibodyModel(
        bodies=[mass],
        joints=[j_gs],
        forces=[f_sd],
    )


if __name__ == "__main__":
    from pmd.gui import PostProcessor

    model = build_model()
    T, uT = model.solve(
        method='CASADI-DAE',
        t_final=T_FINAL,
        t_eval=np.linspace(0, T_FINAL, N_EVAL),
        ic_correct=IC_CORRECT,
    )
    PostProcessor(model=model, T=T, uT=uT).show()
