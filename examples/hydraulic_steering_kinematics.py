"""HSK (Hydraulic Steering Kinematics) - Planar Multibody Dynamics (pmd)"""

import numpy as np

from pmd.core import (
    Actuator,
    Body,
    Function,
    Ground,
    PlanarMultibodyModel,
    RevJoint,
    TranJoint,
    Weight,
)


T_FINAL    = 2.0
N_EVAL     = 2001
IC_CORRECT = True
DELTA_S    = 0.020

P_yoke_pivot = np.array([ 0.00,  0.00])
P_yoke_rod   = np.array([-0.10, -0.28])
P_barrel_gnd = np.array([ 0.40, -0.20])
P_wheel      = np.array([-0.35, -0.05])

L_barrel = 0.30
L_rod    = 0.30

# Inner endpoints of barrel/rod along the cylinder axis (P_barrel_gnd → P_yoke_rod)
_axis          = P_yoke_rod - P_barrel_gnd
L0_cyl         = float(np.linalg.norm(_axis))
_u             = _axis / L0_cyl
P_barrel_inner = P_barrel_gnd + L_barrel * _u
P_rod_inner    = P_yoke_rod   - L_rod    * _u

yoke_wheel, mk_y_pivot, _, mk_y_rod = Body.from_plate(
    P_yoke_pivot, P_wheel, P_yoke_rod,
    mass=31.31, inertia=5.28E+05, name='yoke_wheel')

barrel, mk_b_ground, mk_b_tran = Body.from_link(
    P_barrel_gnd, P_barrel_inner,
    mass=5.0, inertia=0.05, thickness=0.030,
    marker_theta=0.0, name='barrel')

rod, mk_r_tran, mk_r_yoke = Body.from_link(
    P_rod_inner, P_yoke_rod,
    mass=2.0, inertia=0.01, thickness=0.015,
    marker_theta=0.0, name='rod')

g_yoke_pivot = Ground.add_marker(P_yoke_pivot.tolist())
g_barrel_pin = Ground.add_marker(P_barrel_gnd.tolist())

j_gy = RevJoint(iMarker=g_yoke_pivot, jMarker=mk_y_pivot,  name='ground_yoke')
j_gb = RevJoint(iMarker=g_barrel_pin, jMarker=mk_b_ground, name='ground_barrel')
j_yr = RevJoint(iMarker=mk_y_rod,     jMarker=mk_r_yoke,   name='yoke_rod')
j_br = TranJoint(iMarker=mk_b_tran,   jMarker=mk_r_tran,   name='barrel_rod')

S_law = Function(type='b', t_start=0.0, t_end=T_FINAL,
                 f_start=L0_cyl, f_end=L0_cyl + DELTA_S)

act = Actuator(iMarker=mk_r_yoke, jMarker=g_barrel_pin,
               control='length', law=S_law, name='hydraulic_actuator')

fw = Weight(gravity=9.80665)

model = PlanarMultibodyModel(
    bodies=[yoke_wheel, barrel, rod],
    joints=[j_gy, j_gb, j_yr, j_br],
    forces=[fw, act],
    functions=[S_law],
)


if __name__ == "__main__":
    from pmd.gui import PostProcessor

    T, uT = model.solve(
        analysis='kinematic',
        t_final=T_FINAL,
        t_eval=np.linspace(0, T_FINAL, N_EVAL),
        ic_correct=IC_CORRECT,
    )

    PostProcessor(model=model, T=T, uT=uT).show()