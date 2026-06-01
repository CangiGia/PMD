"""
CS (Crank-Slider) Model - Planar Multibody Dynamics (pmd)
=========================================================
KINEMATIC SCHEME & INITIAL STATE
--------------------------------

               (O2)
               *     *
              *            *
             *                   *
            *                          *
           *                                 *     +--------------+
        (O1) - - - - - - - - - - - - - - - - - (O3)|              | - - -> X (Axis)
                                                   +--------------+
                                                   
BODY PARAMETERS
---------------                                       
+------------+------------------+----------------------+-------------------+
| Body Name  | Kinematic Role   | Geometry / Length    |  Initial Angle    |
|            |                  |        [m]           |       [deg]       |
+------------+------------------+----------------------+-------------------+
| crankshaft |      Crank       |     L = 0.10         |      45.0000      |
| rod        |  Connecting Rod  |     L = 0.40         |     −10.1821      |
| slider     |   Slider Block   | Sphere (R = 0.02)*   |       0.0000      |
+------------+------------------+----------------------+-------------------+
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

crankshaft = Body(
    mass=0.8169313258,
    inertia=1.294122161962e-4,
    position=[0.03535533906, 0.03535533906],
    orientation=np.deg2rad(45.0),
    name="crankshaft",
    shape=Link(length=0.10, thickness=0.015),
)

rod = Body(
    mass=2.689171326,
    inertia=4.414522158e-4,
    position=[0.2675608549, 0.03535522707],
    orientation=np.deg2rad(-10.1821),
    name="rod",
    shape=Link(length=0.40, thickness=0.015),
)

slider = Body(
    mass=0.2614140191,
    inertia=4.182624305e-5,
    position=[0.4644110316, 0.0],
    name="slider",
    shape=Rectangle(width=0.08, height=0.04),
)

mk_g_crank = Ground.add_marker([0.0, 0.0])  # crank pivot O1
mk_g_slider = Ground.add_marker([0.4644110316, 0.0], theta=0.0)  # slider guide

mk_c_ground = crankshaft.add_marker([-0.05, 0.0])
mk_c_rod = crankshaft.add_marker([0.05, 0.0])

mk_r_crank = rod.add_marker([-0.2, 0.0])
mk_r_slider = rod.add_marker([0.2, 0.0])

mk_s_rod = slider.add_marker([0.0, 0.0])
mk_s_ground = slider.add_marker([0.0, 0.0], theta=0.0)

fn_mot = Function(
    type="a",
    coeff=[
        np.deg2rad(55.1821),
        np.deg2rad(20.0),
        0.0,
    ],
)

j_gc = RevJoint(iMarker=mk_g_crank, jMarker=mk_c_ground, name="ground_crankshaft_joint")
j_cr = RevJoint(iMarker=mk_c_rod, jMarker=mk_r_crank, name="crankshaft_rod_joint")
j_rs = RevJoint(iMarker=mk_s_rod, jMarker=mk_r_slider, name="rod_slider_joint")
j_gs = TranJoint(iMarker=mk_s_ground, jMarker=mk_g_slider, name="ground_slider_joint")
j_mot = RotMotion(iBody=crankshaft, jBody=rod, iFunct=fn_mot, name="crank_rod_motor")

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
