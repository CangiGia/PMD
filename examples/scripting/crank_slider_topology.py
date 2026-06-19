"""
CS (Crank-Slider) Model — Topology-Driven Assembly
====================================================

KINEMATIC SCHEME & INITIAL STATE
---------------------------------

               (O2)
               *     *
              *            *
             *                   *
            *                          *
           *                                 *     +--------------+
        (O1) - - - - - - - - - - - - - - - - - (O3)|              | - - -> X (Axis)
                                                   +--------------+

JOINT INITIAL VALUES (q0)
--------------------------
+-----------------------------+------------+----------------------------------+
| Joint                       | Type       | q0                               |
+-----------------------------+------------+----------------------------------+
| ground_crankshaft_joint     | RevJoint   | +45.0000 deg  (absolute)         |
| crankshaft_rod_joint        | RevJoint   | −55.1821 deg  (relative to crank)|
| ground_slider_joint         | TranJoint  |   0.4644 m    (slide distance)   |
| rod_slider_joint            | RevJoint   | loop-closing  — no q0 required   |
+-----------------------------+------------+----------------------------------+
"""

import numpy as np

from pmd.core import (
    Body,
    Function,
    Ground,
    PlanarMultibodyModel,
    RevJoint,
    RotMotion,
    TranJoint,
    Weight,
)
from pmd.core.shapes import Link, Rectangle

T_FINAL = 20.0
N_EVAL = 20000
IC_CORRECT = True

# bodies
crankshaft = Body(
    mass=0.8169313258,
    inertia=1.294122161962e-4,
    name="crankshaft",
    shape=Link(length=0.10, thickness=0.015),
)

rod = Body(
    mass=2.689171326,
    inertia=4.414522158e-4,
    name="rod",
    shape=Link(length=0.40, thickness=0.015),
)

slider = Body(
    mass=0.2614140191,
    inertia=4.182624305e-5,
    name="slider",
    shape=Rectangle(width=0.08, height=0.04),
)

# markers
mk_g_crank  = Ground.add_marker([0.0, 0.0])            # crank pivot O1
mk_g_slider = Ground.add_marker([0.0, 0.0], theta=0.0) # origin of the horizontal slide track

mk_c_ground = crankshaft.add_marker([-0.05, 0.0])
mk_c_rod    = crankshaft.add_marker([ 0.05, 0.0])

mk_r_crank  = rod.add_marker([-0.2, 0.0])
mk_r_slider = rod.add_marker([ 0.2, 0.0])

mk_s_rod    = slider.add_marker([0.0, 0.0])
mk_s_ground = slider.add_marker([0.0, 0.0], theta=0.0)

# motion
fn_mot = Function(
    type="a",
    coeff=[
        np.deg2rad(55.1821),   # initial relative angle crankshaft → rod [rad]
        np.deg2rad(20.0),      # angular velocity [rad/s]
        0.0,
    ],
)

# joints
j_gc = RevJoint(
    iMarker=mk_g_crank,
    jMarker=mk_c_ground,
    q0=np.deg2rad(45.0),       # absolute orientation of crankshaft at t = 0
    name="ground_crankshaft_joint",
)

j_cr = RevJoint(
    iMarker=mk_c_rod,
    jMarker=mk_r_crank,
    q0=np.deg2rad(-55.1821),   # rod orientation relative to crankshaft at t = 0
    name="crankshaft_rod_joint",
)

j_rs = RevJoint(
    iMarker=mk_s_rod,
    jMarker=mk_r_slider,
    name="rod_slider_joint",   # loop-closing joint — assembler closes it via NR
)

j_gs = TranJoint(
    iMarker=mk_s_ground,
    jMarker=mk_g_slider,
    q0=0.4644110316,           # initial slide distance from track origin [m]
    name="ground_slider_joint",
)

j_mot = RotMotion(
    iBody=crankshaft,
    jBody=rod,
    iFunct=fn_mot,
    name="crank_rod_motor",
)

fw = Weight(gravity=9.80665)

model = PlanarMultibodyModel(
    bodies=[crankshaft, rod, slider],
    joints=[j_gc, j_cr, j_rs, j_gs, j_mot],
    forces=[fw],
    functions=[fn_mot],
)


if __name__ == "__main__":
    from pmd.gui import PostProcessor

    T, uT = model.solve(
        t_final=T_FINAL,
        t_eval=np.linspace(0, T_FINAL, N_EVAL),
        ic_correct=IC_CORRECT,
    )
    PostProcessor(model=model, T=T, uT=uT).show()
