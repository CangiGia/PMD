"""
FBL (Four-Bar Linkage) Model - Planar Multibody Dynamics (pmd)
==============================================================
KINEMATIC SCHEME
----------------
                                                (O3)
                                            *    *
                                        *       *
                                    *          *
                                *             *
                             *               *
                         *                  *
                     *                     *
                 *                        *
             (O2)                        *
             *                          *
            *                          *
           *                          *
          *                          *
       (O1) - - - - - - - - - - -  (O4)

BODY PARAMETERS
----------------
All links are modeled as standard rigid bodies.

+----------+-----------------+----------+----------+---------------+
| Body     | Kinematic Role  |  Length  |  Width   | Initial Angle |
|          |                 |   [mm]   |  [mm]    |    [deg]      |
+----------+-----------------+----------+----------+---------------+
| link_1   |      Crank      |   40.0   |   8.0    |    40.000     |
| link_2   |     Coupler     |  120.0   |   8.0    |    20.298     |
| link_3   |    Follower     |   80.0   |   8.0    |  −122.605     |
+----------+-----------------+----------+----------+---------------+

TOPOLOGY & JOINTS
-----------------
* O1 [link_1_ground_joint] : RevJoint     Ground <-> link_1  @ (0, 0)
* O2 [link_1_link_2_joint] : RevJoint     link_1 <-> link_2  @ Crank tip
* O3 [j_23]                : RevJoint     link_2 <-> link_3  @ Coupler tip
* O4 [j_3g]                : RevJoint     link_3 <-> Ground  @ (~100, 0) [Loop Closure]

ACTIVE FORCES
-------------
* Driver (motion) : RotMotion imposing 25 deg/s constant angular velocity on 
                    link_1_ground_joint.
* Gravity (fw)   : Active along the -Y axis.
"""

import numpy as np
from pmd.core import (
    Body,
    Function,
    Ground,
    PlanarMultibodyModel,
    RotMotion,
    RevJoint,
    Weight,
)
from pmd.core.shapes import Link
from pmd.core.model import _GroundType


T_FINAL = 25.0
N_EVAL = 25001
IC_CORRECT = True

# Geometry
_L1 = 40.0e-3    # m  crank length
_L2 = 120.0e-3   # m  coupler length
_L3 = 80.0e-3    # m  follower length

_THETA1_DEG = 40.000
_THETA2_DEG = 20.298
_THETA3_DEG = -122.605

_O4 = [0.100, 0.000]
_OMEGA_DEG = -25.0

if _GroundType._instance is not None:
    _GroundType._markers = [_GroundType._instance.origin]

link_1 = Body(mass=1.1528090607e-02, inertia=7.4696175803e-05, name='link_1',
              shape=Link(length=_L1, thickness=8.0e-3))
link_2 = Body(mass=3.1498650607e-02, inertia=0.2078332421e-05, name='link_2',
              shape=Link(length=_L2, thickness=8.0e-3))
link_3 = Body(mass=2.1513370607e-02, inertia=0.141264709e-03,  name='link_3',
              shape=Link(length=_L3, thickness=8.0e-3))

g_l1 = Ground.add_marker([0.000, 0.0])  # crank pivot O1
g_l3 = Ground.add_marker(_O4)           # follower pivot O4

l1_g  = link_1.add_marker([-_L1 / 2, 0.0])
l1_l2 = link_1.add_marker([ _L1 / 2, 0.0])

l2_l1 = link_2.add_marker([-_L2 / 2, 0.0])
l2_l3 = link_2.add_marker([ _L2 / 2, 0.0])

l3_l2 = link_3.add_marker([-_L3 / 2, 0.0])
l3_g  = link_3.add_marker([ _L3 / 2, 0.0])

fn_mot = Function(
    type='a',
    coeff=[
        np.deg2rad(_THETA1_DEG),
        np.deg2rad(_OMEGA_DEG),
        0.0,
    ],
)

j_l1_g = RevJoint(iMarker=g_l1,  jMarker=l1_g,
                  q0=np.deg2rad(_THETA1_DEG),
                  name='link_1_ground_joint')
j_l1_l2 = RevJoint(iMarker=l1_l2, jMarker=l2_l1,
                   q0=np.deg2rad(_THETA2_DEG - _THETA1_DEG),
                   name='link_1_link_2_joint')
j_l2_l3 = RevJoint(iMarker=l2_l3, jMarker=l3_l2,
                   q0=np.deg2rad(_THETA3_DEG - _THETA2_DEG),
                   name='link_2_link_3_joint')
j_l3_g = RevJoint(iMarker=l3_g,  jMarker=g_l3,
                  q0=np.deg2rad(_THETA3_DEG),
                  name='link_3_ground_joint')
motion = RotMotion(iBody=link_1, jBody=Ground, iFunct=fn_mot,
                     name='crank_motion')

fw = Weight(gravity=9.80665)

model = PlanarMultibodyModel(
    bodies=[link_1, link_2, link_3],
    joints=[j_l1_g, j_l1_l2, j_l2_l3, j_l3_g, motion],
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
