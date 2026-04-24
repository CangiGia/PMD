"""
FBL (Four-Bar Linkage) Model - Planar Multibody Dynamics (PMD)
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
* Driver (motion) : RelRotJoint imposing 25 deg/s constant angular velocity on 
                    link_1_ground_joint.
* Gravity (fw)   : Active along the -Y axis.
"""

import numpy as np
from PMD.src import *
from gui.app import PostProcessor


L1 = 40.0e-3    # m  crank length
L2 = 120.0e-3   # m  coupler length
L3 = 80.0e-3    # m  follower length

THETA1_DEG = 40.000    # link_1 IC angle [deg]
THETA2_DEG = 20.298    # link_2 IC angle [deg]
THETA3_DEG = -122.605  # link_3 IC angle [deg]

O4 = [0.100, 0.000]  # m
OMEGA_DEG = -25.0    # deg/s


#%% Bodies
link_1 = Body(
    mass=1.1528090607e-02,     # kg
    inertia=7.4696175803e-05,  # kg*m^2
    name='link_1',
)

link_2 = Body(
    mass=3.1498650607e-02,     # kg
    inertia=0.2078332421e-05,  # kg*m^2
    name='link_2',
)

link_3 = Body(
    mass=2.1513370607e-02,     # kg
    inertia=0.141264709e-03,   # kg*m^2
    name='link_3',
)


#%% Markers
# Ground markers
ground_link_1_ref_frame = Ground.add_marker([0.000, 0.0])  # crank pivot   O1 = (0, 0)
ground_link_3_ref_frame = Ground.add_marker(O4)            # follower pivot O4 ≈ (100 mm, 0)

# link_1 markers (half-length from CM along local x-axis)
link_1_ground_ref_frame = link_1.add_marker([-L1 / 2, 0.0])  # ground-side end
link_1_link_2_ref_frame = link_1.add_marker([L1 / 2, 0.0])   # coupler-side end

# link_2 markers (half-length from CM along local x-axis)
link_2_link_1_ref_frame = link_2.add_marker([-L2 / 2, 0.0])  # link_1-side end
link_2_link_3_ref_frame = link_2.add_marker([L2 / 2, 0.0])   # link_3-side end

# link_3 markers (half-length from CM along local x-axis)
link_3_link_2_ref_frame = link_3.add_marker([-L3 / 2, 0.0])  # link_2-side end
link_3_ground_ref_frame = link_3.add_marker([L3 / 2, 0.0])   # ground-side end


#%% Function (motion driver)
# phi_link1 = c0 + c1*t
fn_mot = Function(
    type='a',
    coeff=[
        np.deg2rad(THETA1_DEG),  # c0: initial absolute angle [rad]
        np.deg2rad(OMEGA_DEG),   # c1: angular velocity [rad/s]
        0.0,                     # c2
    ],
)


#%% Joints
link_1_ground_joint = RevJoint(
    iMarker=ground_link_1_ref_frame,
    jMarker=link_1_ground_ref_frame,
    q0=np.deg2rad(THETA1_DEG),
    name='link_1_ground_joint',
)

link_1_link_2_joint = RevJoint(
    iMarker=link_1_link_2_ref_frame,
    jMarker=link_2_link_1_ref_frame,
    q0=np.deg2rad(THETA2_DEG - THETA1_DEG),
    name='link_1_link_2_joint',
)

link_2_link_3_joint = RevJoint(
    iMarker=link_2_link_3_ref_frame,
    jMarker=link_3_link_2_ref_frame,
    q0=np.deg2rad(THETA3_DEG - THETA2_DEG),
    name='link_2_link_3_joint',
)

link_3_ground_joint = RevJoint(  # loop-closing
    iMarker=link_3_ground_ref_frame,
    jMarker=ground_link_3_ref_frame,
    q0=np.deg2rad(THETA3_DEG),
    name='link_3_ground_joint',
)

motion = RelRotJoint(
    iBody=link_1,
    jBody=Ground,
    iFunct=fn_mot,
    name='crank_motion',
)


#%% Forces
fw = Weight(gravity=9.80665)  # m/s^2, along -y


#%% Simulation
model = PlanarMultibodyModel(
    bodies=[link_1, link_2, link_3],
    joints=[
        link_1_ground_joint,
        link_1_link_2_joint,
        link_2_link_3_joint,
        link_3_ground_joint,
        motion,
    ],
    forces=[fw],
    functions=[fn_mot],
)

T, uT = model.solve(
    method='CASADI-DAE',
    t_final=25,
    t_eval=np.linspace(0, 25, 25001),
    ic_correct=True,
)

post_proc = PostProcessor(model=model, T=T, uT=uT)
post_proc.show()