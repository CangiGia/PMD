"""MSD (Mass-Spring-Damper) Model - Planar Multibody Dynamics (pmd)

          +---------+
          |         |
          |    m    |  mass  (CM design at y = 1.05 m, initial at y = 1.25 m)
          |         |
          +---------+  <- mk_m_spring (y_local = -0.25 m)
            |     |
            /     |
            \    _|_
    spring  /   |___|
            \   | | | damper
            /     |
            |     |
         ===+=====+=== ground (y = 0.0 m)

ADAMS REFERENCE PARAMETERS (converted from mm, N, kg, s to SI: m, N, kg, s)
---------------------------------------------------------------------------
+------+-----------+---------------+-----------------+--------------------+
| Body |  m  [kg]  |  I  [kg·m²]   | Design CM [m]   | Initial CM (t = 0) |
+------+-----------+---------------+-----------------+--------------------+
| mass |  10000.0  |  1.0e-3       | (0.0, 1.05)     | (0.0, 1.25)        |
+------+-----------+---------------+-----------------+--------------------+

JOINT & FORCE PARAMETERS
------------------------
- Translational joint along +Y at design CM (0.0, 1.05 m) with initial
  displacement q0 = +0.20 m (200 mm in ADAMS -> initial spring stretch)
- Spring-damper between ground (0.0, 0.0 m) and mass bottom face
  (local [0.0, -0.25 m]):
    k  = 1.0e6 N/m       (1.0e3 N/mm)
    dc = 1.0e4 N·s/m     (10.0 N·s/mm)
    L0 = 0.80 m          (800 mm design distance, zero preload at design pose)
- Gravity: disabled
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

T_FINAL = 15.0
N_EVAL = 3001
IC_CORRECT = True

# Geometry & physical parameters (SI units: m, kg, N, s)
_BLOCK_SIDE = 0.50       # 500 mm block side
_Y_CM_DESIGN = 1.05      # 1050 mm design CM height (800 mm corner + 250 mm half-side)
_Q0_TRAN = 0.20          # 200 mm initial translational joint displacement (translational_ic)
_L0_SPRING = 0.80        # 800 mm undeformed spring length (design distance DM)
_K_SPRING = 1.0e6        # 1.0e3 N/mm = 1.0e6 N/m
_C_DAMPER = 1.0e4        # 10 N·s/mm  = 1.0e4 N·s/m

# Body (no explicit position/orientation -> assembled from TranJoint q0)
mass = Body(
    mass=10000.0,
    inertia=1.0e-3,
    name="mass",
    shape=Rectangle(width=_BLOCK_SIDE, height=_BLOCK_SIDE),
)

# Ground markers
mk_g_spring = Ground.add_marker([0.0, 0.0], name="ground_spring_marker")
mk_g_tran = Ground.add_marker(
    [0.0, _Y_CM_DESIGN],
    theta=np.pi / 2,
    name="ground_tran_marker",
)

# Mass markers (local coordinates relative to mass CM)
mk_m_spring = mass.add_marker(
    [0.0, -_BLOCK_SIDE / 2],
    name="mass_spring_marker",
)
mk_m_tran = mass.add_marker(
    [0.0, 0.0],
    theta=np.pi / 2,
    name="mass_tran_marker",
)

# Translational joint with initial displacement q0 = 0.20 m
# Places mass CM at y(0) = _Y_CM_DESIGN + _Q0_TRAN = 1.25 m
j_gs = TranJoint(
    iMarker=mk_m_tran,
    jMarker=mk_g_tran,
    q0=_Q0_TRAN,
    name="mass_ground_joint",
)

# Point-to-point spring-damper (stretched to L(0) = 1.00 m at t = 0, L0 = 0.80 m)
f_sd = PtpForce(
    iMarker=mk_m_spring,
    jMarker=mk_g_spring,
    k=_K_SPRING,
    L0=_L0_SPRING,
    dc=_C_DAMPER,
    name="spring",
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
        t_eval=np.linspace(0.0, T_FINAL, N_EVAL),
        ic_correct=IC_CORRECT,
    )
    PostProcessor(model=model, T=T, uT=uT).show()

