"""HSK (Hydraulic Steering Kinematics) - Planar Multibody Dynamics (pmd)"""

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
L_CYL = 797.306  # hydraulic full extended length
L0_CYL = 680.646  # hydraulic cylinder initial length
DELTA_S = L_CYL - L0_CYL  # hydraulic cylinder available stroke

# bodies
# barrel, barrel_ground_reference_frame, barrel_rod_reference_frame = Body.as_link(
#     np.array([0.0, 0.0]), np.array([-439.957,-45.106]),
#     mass=4.55, inertia=86113.72, thickness=77,
#     marker_theta=0.0, name='barrel')

# rod, rod_end_reference_frame, rod_yoke_reference_frame = Body.as_link(
#     np.array([-212.893, -33.092]), np.array([-676.982, -70.529]),
#     mass=4.52, inertia=81024.95, thickness=40,
#     marker_theta=0.0, name='rod')

barrel, barrel_ground_reference_frame, barrel_rod_reference_frame = Body.as_link(
    np.array([0.0, 0.0]),
    np.array([-432.957, -45.106]),
    mass=4.55,
    inertia=86113.72,
    thickness=77,
    marker_theta=0.0,
    name="barrel",
)

rod, rod_end_reference_frame, rod_yoke_reference_frame = Body.as_link(
    np.array([-207.733, -21.642]),
    np.array([-676.982, -70.529]),
    mass=4.52,
    inertia=81024.95,
    thickness=40,
    marker_theta=0.0,
    name="rod",
)

yoke, yoke_rod_reference_frame, yoke_ground_reference_frame, yoke_wheel_reference_frame = (
    Body.as_plate(
        np.array([-676.982, -70.529]),
        np.array([-712.987, 215.573]),
        np.array([-905.3226, 215.573]),
        mass=31.31,
        inertia=5.28e05,
        name="yoke",
    )
)

wheel_cm = np.array([-905.3226, 215.573])
wheel = Body(
    position=wheel_cm,
    orientation=0.0,
    shape=Rectangle(width=300.0, height=412.0),
    mass=127,
    inertia=2.695e06,
    name="wheel",
)

# markers
ground_yoke_reference_frame = Ground.add_marker(np.array([-712.987, 215.573]))
ground_barrel_reference_frame = Ground.add_marker(np.array([0.0, 0.0]))

# joints
ground_yoke_revolute_joint = RevJoint(
    iMarker=ground_yoke_reference_frame,
    jMarker=yoke_ground_reference_frame,
    name="ground_yoke_rev_joint",
)

ground_barrel_revolute_joint = RevJoint(
    iMarker=ground_barrel_reference_frame,
    jMarker=barrel_ground_reference_frame,
    name="ground_barrel_rev_joint",
)

rod_yoke_revolute_joint = RevJoint(
    iMarker=rod_yoke_reference_frame, jMarker=yoke_rod_reference_frame, name="yoke_rod_rev_joint"
)

barrel_rod_translational_joint = TranJoint(
    iMarker=barrel_rod_reference_frame,
    jMarker=rod_end_reference_frame,
    name="barrel_rod_tran_joint",
)

yoke_wheel_rigid_joint = RevJoint(
    iMarker=yoke_wheel_reference_frame, jMarker=wheel.cm, fix=1, name="yoke_wheel_rigid_joint"
)

# forces and actuators
extension_law = Function(
    type="b",
    t_start=0.0,
    t_end=T_FINAL,
    f_start=L0_CYL,
    f_end=L0_CYL + DELTA_S,
)

actuator = Actuator(
    iMarker=rod_yoke_reference_frame,
    jMarker=ground_barrel_reference_frame,
    control="length",
    law=extension_law,
)

model = PlanarMultibodyModel(
    bodies=[yoke, barrel, rod, wheel],
    joints=[
        ground_yoke_revolute_joint,
        ground_barrel_revolute_joint,
        rod_yoke_revolute_joint,
        barrel_rod_translational_joint,
        yoke_wheel_rigid_joint,
    ],
    forces=[actuator],
    functions=[extension_law],
    units=units,
)

if __name__ == "__main__":
    from pmd.gui import PostProcessor

    T, uT = model.solve(
        analysis="kinematic",
        t_final=T_FINAL,
        t_eval=np.linspace(0, T_FINAL, N_EVAL),
        ic_correct=IC_CORRECT,
    )

    PostProcessor(model=model, T=T, uT=uT).show()
