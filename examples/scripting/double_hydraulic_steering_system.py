"""
Double HSS (Hydraulic Steering System) - Kinematic analysis - 
Planar Multibody Dynamics (pmd)
"""

import numpy as np

from pmd.core import (
    Actuator,
    Body,
    Function,
    Ground,
    PlanarMultibodyModel,
    Rectangle,
    RevJoint,
    TranJoint,
    UnitSystem,
)

units = UnitSystem(length="mm", force="N", angle="rad")

T_FINAL = 2.0  # assumption, all the flow used for the cylinder extension
N_EVAL = int(T_FINAL / 0.001) + 1
IC_CORRECT = True
L_CYL = 787.20   # hydraulic full extended length (kinematic limit of updated geometry)
L0_CYL = 680.646  # hydraulic cylinder initial length
DELTA_S = L_CYL - L0_CYL  # hydraulic cylinder available stroke
SHIFT = 130.026  # x-offset of the right-side tie-rod yoke attachment

# =============================================================================
# Tie Rod
# =============================================================================

tie_rod, tie_rod_yoke_left_reference_frame, tie_rod_yoke_right_reference_frame = Body.as_link(
    np.array([-610.494, -103.872]),
    np.array([610.494 + SHIFT, -103.872]),
    mass=6.74,
    inertia=1234456.10,
    thickness=35,
    marker_theta=0.0,
    name="tie_rod",
    color="tab:red",
)

# =============================================================================
# LEFT side  (X < 0)
# =============================================================================

# bodies
barrel_left, barrel_ground_reference_frame_left, barrel_rod_reference_frame_left = Body.as_link(
    np.array([0.0, 0.0]),
    np.array([-423.011, -44.070]),
    mass=4.55,
    inertia=86113.72,
    thickness=77,
    marker_theta=0.0,
    name="barrel_left",
    color="tab:orange",
)

rod_left, rod_end_reference_frame_left, rod_yoke_reference_frame_left = Body.as_link(
    np.array([-207.733, -21.642]),
    np.array([-676.982, -70.529]),
    mass=4.52,
    inertia=81024.95,
    thickness=40,
    marker_theta=0.0,
    name="rod_left",
    color="tab:gray",
)

(
    yoke_left,
    yoke_rod_reference_frame_left,
    yoke_ground_reference_frame_left,
    yoke_wheel_reference_frame_left,
) = Body.as_plate(
    np.array([-676.982, -70.529]),
    np.array([-712.987, 142.348]),
    np.array([-905.3226, 142.348]),
    mass=31.31,
    inertia=5.28e05,
    name="yoke_left",
    color="tab:blue",
)

wheel_left = Body(
    position=np.array([-905.3226, 142.348]),
    orientation=0.0,
    shape=Rectangle(width=300.0, height=412.0),
    mass=127,
    inertia=2.695e06,
    name="wheel_left",
    color="silver",
)

# ground markers
ground_yoke_reference_frame_left = Ground.add_marker(np.array([-712.987, 142.348]))
ground_barrel_reference_frame_left = Ground.add_marker(np.array([0.0, 0.0]))

# deferred marker: tie-rod attachment on yoke_left
yoke_left_tie_rod_reference_frame = yoke_left.add_marker_at(
    tie_rod_yoke_left_reference_frame, theta=0.0
)

# joints
ground_yoke_revolute_joint_left = RevJoint(
    iMarker=ground_yoke_reference_frame_left,
    jMarker=yoke_ground_reference_frame_left,
    name="ground_yoke_rev_joint_left",
)

ground_barrel_revolute_joint_left = RevJoint(
    iMarker=ground_barrel_reference_frame_left,
    jMarker=barrel_ground_reference_frame_left,
    name="ground_barrel_rev_joint_left",
)

rod_yoke_revolute_joint_left = RevJoint(
    iMarker=rod_yoke_reference_frame_left,
    jMarker=yoke_rod_reference_frame_left,
    name="yoke_rod_rev_joint_left",
)

barrel_rod_translational_joint_left = TranJoint(
    iMarker=barrel_rod_reference_frame_left,
    jMarker=rod_end_reference_frame_left,
    name="barrel_rod_tran_joint_left",
)

yoke_wheel_rigid_joint_left = RevJoint(
    iMarker=yoke_wheel_reference_frame_left,
    jMarker=wheel_left.cm,
    fix=1,
    name="yoke_wheel_rigid_joint_left",
)

# =============================================================================
# RIGHT side  (X > 0)
# =============================================================================

# bodies
barrel_right, barrel_ground_reference_frame_right, barrel_rod_reference_frame_right = Body.as_link(
    np.array([SHIFT, 0.0]),
    np.array([423.011 + SHIFT, -44.070]),
    mass=4.55,
    inertia=86113.72,
    thickness=77,
    marker_theta=0.0,
    name="barrel_right",
    color="tab:orange",
)

rod_right, rod_end_reference_frame_right, rod_yoke_reference_frame_right = Body.as_link(
    np.array([207.733 + SHIFT, -21.642]),
    np.array([676.982 + SHIFT, -70.529]),
    mass=4.52,
    inertia=81024.95,
    thickness=40,
    marker_theta=0.0,
    name="rod_right",
    color="tab:gray",
)

(
    yoke_right,
    yoke_wheel_reference_frame_right,
    yoke_ground_reference_frame_right,
    yoke_rod_reference_frame_right,
) = Body.as_plate(
    np.array([905.3226 + SHIFT, 142.348]),
    np.array([712.987 + SHIFT, 142.348]),
    np.array([676.982 + SHIFT, -70.529]),
    mass=31.31,
    inertia=5.28e05,
    name="yoke_right",
    color="tab:blue",
)

wheel_right = Body(
    position=np.array([905.3226 + SHIFT, 142.348]),
    orientation=0.0,
    shape=Rectangle(width=300.0, height=412.0),
    mass=127,
    inertia=2.695e06,
    name="wheel_right",
    color="silver",
)

# ground markers
ground_yoke_reference_frame_right = Ground.add_marker(np.array([712.987 + SHIFT, 142.348]))
ground_barrel_reference_frame_right = Ground.add_marker(np.array([SHIFT, 0.0]))

# deferred marker: tie-rod attachment on yoke_right
yoke_right_tie_rod_reference_frame = yoke_right.add_marker_at(
    tie_rod_yoke_right_reference_frame, theta=0.0
)

# joints
ground_yoke_revolute_joint_right = RevJoint(
    iMarker=ground_yoke_reference_frame_right,
    jMarker=yoke_ground_reference_frame_right,
    name="ground_yoke_rev_joint_right",
)

ground_barrel_revolute_joint_right = RevJoint(
    iMarker=ground_barrel_reference_frame_right,
    jMarker=barrel_ground_reference_frame_right,
    name="ground_barrel_rev_joint_right",
)

rod_yoke_revolute_joint_right = RevJoint(
    iMarker=rod_yoke_reference_frame_right,
    jMarker=yoke_rod_reference_frame_right,
    name="yoke_rod_rev_joint_right",
)

barrel_rod_translational_joint_right = TranJoint(
    iMarker=barrel_rod_reference_frame_right,
    jMarker=rod_end_reference_frame_right,
    name="barrel_rod_tran_joint_right",
)

yoke_wheel_rigid_joint_right = RevJoint(
    iMarker=yoke_wheel_reference_frame_right,
    jMarker=wheel_right.cm,
    fix=1,
    name="yoke_wheel_rigid_joint_right",
)

# =============================================================================
# Tie-rod revolute joints (connect tie_rod to both yokes)
# =============================================================================

yoke_right_tie_rod_revolute_joint = RevJoint(
    iMarker=yoke_right_tie_rod_reference_frame,
    jMarker=tie_rod_yoke_right_reference_frame,
    name="yoke_right_tie_rod_rev_joint",
)

yoke_left_tie_rod_revolute_joint = RevJoint(
    iMarker=yoke_left_tie_rod_reference_frame,
    jMarker=tie_rod_yoke_left_reference_frame,
    name="yoke_left_tie_rod_rev_joint",
)

# =============================================================================
# Forces and actuators  —  both sides extend simultaneously
# =============================================================================

extension_law_left = Function(
    type="b",
    t_start=0.0,
    t_end=T_FINAL,
    f_start=L0_CYL,
    f_end=L0_CYL + DELTA_S,
)

extension_law_right = Function(
    type="b",
    t_start=0.0,
    t_end=T_FINAL,
    f_start=L0_CYL,
    f_end=L0_CYL + DELTA_S,
)

actuator_left = Actuator(
    iMarker=rod_yoke_reference_frame_left,
    jMarker=ground_barrel_reference_frame_left,
    control="length",
    law=extension_law_left,
)

actuator_right = Actuator(
    iMarker=rod_yoke_reference_frame_right,
    jMarker=ground_barrel_reference_frame_right,
    control="length",
    law=extension_law_right,
)

# =============================================================================
# Model
# =============================================================================

double_hydralic_steering_model = PlanarMultibodyModel(
    bodies=[
        yoke_left,
        barrel_left,
        rod_left,
        wheel_left,
        yoke_right,
        barrel_right,
        rod_right,
        wheel_right,
        tie_rod,
    ],
    joints=[
        ground_yoke_revolute_joint_left,
        ground_barrel_revolute_joint_left,
        rod_yoke_revolute_joint_left,
        barrel_rod_translational_joint_left,
        yoke_wheel_rigid_joint_left,
        ground_yoke_revolute_joint_right,
        ground_barrel_revolute_joint_right,
        rod_yoke_revolute_joint_right,
        barrel_rod_translational_joint_right,
        yoke_wheel_rigid_joint_right,
        yoke_right_tie_rod_revolute_joint,
        yoke_left_tie_rod_revolute_joint,
    ],
    # forces=[actuator_left],
    # functions=[extension_law_left],
    forces=[actuator_right],
    functions=[extension_law_right],
    units=units,
)

if __name__ == "__main__":
    from pmd.gui import PostProcessor, preview_model

    # uncomment below to preview the model before running the analysis
    # preview_model(double_hydralic_steering_model)

    T, uT = double_hydralic_steering_model.solve(
        analysis="kinematic",
        t_final=T_FINAL,
        t_eval=np.linspace(0, T_FINAL, N_EVAL),
        ic_correct=IC_CORRECT,
    )

    PostProcessor(model=double_hydralic_steering_model, T=T, uT=uT).show()
