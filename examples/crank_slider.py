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
All links are modeled as standard rigid bodies.

+------------+------------------+----------------------+-------------------+
| Body Name  | Kinematic Role   | Geometry / Length    |  Initial Angle    |
|            |                  |        [m]           |       [deg]       |
+------------+------------------+----------------------+-------------------+
| crankshaft |      Crank       |     L = 0.10         |      45.0000      |
| rod        |  Connecting Rod  |     L = 0.40         |     −10.1821      |
| slider     |   Slider Block   | Sphere (R = 0.02)*   |       0.0000      |
+------------+------------------+----------------------+-------------------+

TOPOLOGY & JOINTS
-----------------
* O1 [j_gc]  : RevJoint   Ground <-> Crankshaft  @ (0, 0)
* O2 [j_cr]  : RevJoint   Crankshaft <-> Rod     @ Crank tip
* O3 [j_rs]  : RevJoint   Rod <-> Slider         @ Rod tip
* X  [j_gs]  : TranJoint  Ground <-> Slider      @ Y=0 (Sliding along X-axis)

BOUNDARY CONDITIONS & SOLVER
----------------------------
* Driver (j_mot) : RelRotJoint on O2 (j_cr). Imposes a constant relative angular 
                   velocity of 20 deg/s BETWEEN the Crankshaft and the Rod.
* Gravity (fw)   : Active along the -Y axis.
"""

import numpy as np
from pmd.core import (
    Body,
    Function,
    Ground,
    PlanarMultibodyModel,
    RelRotJoint,
    RevJoint,
    TranJoint,
    Weight,
)
from pmd.core.model import _GroundType


T_FINAL = 20.0
N_EVAL = 20000
IC_CORRECT = True


def build_model() -> PlanarMultibodyModel:
    """Build the crank-slider :class:`PlanarMultibodyModel` (no solve)."""
    gt = _GroundType._instance
    if gt is not None:
        _GroundType._markers = [gt.origin]

    # Bodies
    crankshaft = Body(
        mass=0.8169313258,           # kg
        inertia=1.294122161962e-4,   # kg*m^2
        position=[0.03535533906, 0.03535533906],  # m
        orientation=np.deg2rad(45.0),
        name='crankshaft',
    )

    rod = Body(
        mass=2.689171326,            # kg
        inertia=4.414522158e-4,      # kg*m^2
        position=[0.2675608549, 0.03535522707],  # m
        orientation=np.deg2rad(-10.1821),
        name='rod',
    )

    slider = Body(
        mass=0.2614140191,           # kg
        inertia=4.182624305e-5,      # kg*m^2
        position=[0.4644110316, 0.0],  # m
        name='slider',
    )

    # Ground markers
    mk_g_crank = Ground.add_marker([0.0, 0.0])                       # crank pivot O1
    mk_g_slider = Ground.add_marker([0.4644110316, 0.0], theta=0.0)  # slider guide

    # Crankshaft markers (half-length = 0.05 m)
    mk_c_ground = crankshaft.add_marker([-0.05, 0.0])
    mk_c_rod = crankshaft.add_marker([0.05, 0.0])

    # Rod markers (half-length = 0.2 m)
    mk_r_crank = rod.add_marker([-0.2, 0.0])
    mk_r_slider = rod.add_marker([0.2, 0.0])

    # Slider markers
    mk_s_rod = slider.add_marker([0.0, 0.0])
    mk_s_ground = slider.add_marker([0.0, 0.0], theta=0.0)

    # Function: phi_crank - phi_rod = c0 + c1*t
    fn_mot = Function(
        type='a',
        coeff=[
            np.deg2rad(55.1821),
            np.deg2rad(20.0),
            0.0,
        ],
    )

    # Joints
    j_gc = RevJoint(iMarker=mk_g_crank, jMarker=mk_c_ground,
                    name='ground_crankshaft_joint')
    j_cr = RevJoint(iMarker=mk_c_rod, jMarker=mk_r_crank,
                    name='crankshaft_rod_joint')
    j_rs = RevJoint(iMarker=mk_s_rod, jMarker=mk_r_slider,
                    name='rod_slider_joint')
    j_gs = TranJoint(iMarker=mk_s_ground, jMarker=mk_g_slider,
                     name='ground_slider_joint')
    j_mot = RelRotJoint(iBody=crankshaft, jBody=rod, iFunct=fn_mot,
                        name='crank_rod_motor')

    # Forces
    fw = Weight(gravity=9.80665)

    return PlanarMultibodyModel(
        bodies=[crankshaft, rod, slider],
        joints=[j_gc, j_cr, j_rs, j_gs, j_mot],
        forces=[fw],
        functions=[fn_mot],
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
