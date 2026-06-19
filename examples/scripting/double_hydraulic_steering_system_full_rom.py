"""
Double HSS (Hydraulic Steering System) - Full ROM kinematic analysis
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

T_FINAL = 2.0                           # simulation duration for each half-ROM [s]
N_EVAL  = int(T_FINAL / 0.001) + 1     # number of output steps
IC_CORRECT = True
L_CYL   = 787.20                        # fully extended cylinder length [mm]
L0_CYL  = 680.646                       # rest (neutral) cylinder length [mm]
DELTA_S = L_CYL - L0_CYL              # extension stroke [mm]  (retraction is asymmetric — see __main__)
SHIFT   = 130.026                       # x-offset of right-side tie-rod attachment [mm]


# =============================================================================
# Model factory
# =============================================================================

def build_model(f_start: float, f_end: float) -> PlanarMultibodyModel:
    """Construct a complete, independent model instance.

    All objects (bodies, markers, joints, function, actuator) are created
    fresh on every call so that two simulations can be run without any
    shared mutable state.

    Parameters
    ----------
    f_start : float
        Cylinder length imposed at t = 0 [mm].
    f_end : float
        Cylinder length imposed at t = T_FINAL [mm].

    Returns
    -------
    PlanarMultibodyModel
    """

    # ── Tie Rod ──────────────────────────────────────────────────────────────
    (
        tie_rod,
        tie_rod_yoke_left_reference_frame,
        tie_rod_yoke_right_reference_frame,
    ) = Body.as_link(
        np.array([-610.494, -103.872]),
        np.array([610.494 + SHIFT, -103.872]),
        mass=6.74,
        inertia=1234456.10,
        thickness=35,
        marker_theta=0.0,
        name="tie_rod",
        color="tab:red",
    )

    # ── LEFT side (X < 0) ────────────────────────────────────────────────────

    # bodies
    (
        barrel_left,
        barrel_ground_reference_frame_left,
        barrel_rod_reference_frame_left,
    ) = Body.as_link(
        np.array([0.0, 0.0]),
        np.array([-423.011, -44.070]),
        mass=4.55,
        inertia=86113.72,
        thickness=77,
        marker_theta=0.0,
        name="barrel_left",
        color="tab:orange",
    )

    (
        rod_left,
        rod_end_reference_frame_left,
        rod_yoke_reference_frame_left,
    ) = Body.as_link(
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
    ground_yoke_reference_frame_left   = Ground.add_marker(np.array([-712.987, 142.348]))
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

    # ── RIGHT side (X > 0) ───────────────────────────────────────────────────

    # bodies
    (
        barrel_right,
        barrel_ground_reference_frame_right,
        barrel_rod_reference_frame_right,
    ) = Body.as_link(
        np.array([SHIFT, 0.0]),
        np.array([423.011 + SHIFT, -44.070]),
        mass=4.55,
        inertia=86113.72,
        thickness=77,
        marker_theta=0.0,
        name="barrel_right",
        color="tab:orange",
    )

    (
        rod_right,
        rod_end_reference_frame_right,
        rod_yoke_reference_frame_right,
    ) = Body.as_link(
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
    ground_yoke_reference_frame_right   = Ground.add_marker(np.array([712.987 + SHIFT, 142.348]))
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

    # ── Tie-rod revolute joints ───────────────────────────────────────────────
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

    # ── Actuator (right cylinder drives the full system via tie rod) ──────────
    law = Function(
        type="b",
        t_start=0.0,
        t_end=T_FINAL,
        f_start=f_start,
        f_end=f_end,
    )
    actuator = Actuator(
        iMarker=rod_yoke_reference_frame_right,
        jMarker=ground_barrel_reference_frame_right,
        control="length",
        law=law,
    )

    return PlanarMultibodyModel(
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
        forces=[actuator],
        functions=[law],
        units=units,
    )



def _patch_full_rom(model_ret, model_ext, T_full, uT_full) -> None:
    """Overwrite model_ext._result_container with concatenated full-ROM data.

    The PostProcessor reads curve data from each Body/Joint
    _result_container, not from the uT matrix stored in the Session.
    After two separate kinematic solves the containers only hold N_EVAL steps
    each (the last solve wins).  This helper rebuilds them so that every array
    spans the full time history (2*N_EVAL - 1 steps):

      * Positions / velocities  -- sliced from uT_full columns directly.
      * Accelerations           -- concatenated from the two individual RCs.
      * Joint reactions         -- concatenated from the two individual RCs.
    """
    nB  = len(model_ext.Bodies)
    nB3 = 3 * nB

    for Bi, (b_ret, b_ext) in enumerate(zip(model_ret.Bodies, model_ext.Bodies)):
        i      = 3 * Bi
        rc_ret = b_ret._result_container
        rc_ext = b_ext._result_container

        def _cat_acc(key, _r=rc_ret, _e=rc_ext):
            return np.concatenate([
                _r["accelerations"][key][::-1][:-1],
                _e["accelerations"][key],
            ])

        b_ext._result_container = {
            "positions": {
                "x":   uT_full[:, i],
                "y":   uT_full[:, i + 1],
                "phi": uT_full[:, i + 2],
            },
            "velocities": {
                "dx":   uT_full[:, nB3 + i],
                "dy":   uT_full[:, nB3 + i + 1],
                "dphi": uT_full[:, nB3 + i + 2],
            },
            "accelerations": {
                "ddx":   _cat_acc("ddx"),
                "ddy":   _cat_acc("ddy"),
                "ddphi": _cat_acc("ddphi"),
            },
        }

    for j_ret, j_ext in zip(model_ret.Joints, model_ext.Joints):
        rxn_full = np.concatenate([
            j_ret._result_container["reactions"][::-1][:-1],
            j_ext._result_container["reactions"],
        ], axis=0)
        j_ext._result_container = {"reactions": rxn_full}


if __name__ == "__main__":
    from pmd.gui import PostProcessor

    t_eval = np.linspace(0.0, T_FINAL, N_EVAL)

    # -- Simulation B: rest -> maximum extension (steer right) ---------------
    #
    #   Run first: its result is used to compute the true retraction limit
    #   (L_CYL_MIN) for Sim A via the geometric mirror symmetry of the
    #   mechanism (both sides have identical local kinematics).
    #
    print("\n[Full ROM]  Sim B -- extension   (L0_CYL -> L_CYL, steer right)")
    model_ext = build_model(f_start=L0_CYL, f_end=L_CYL)
    T_ext, uT_ext = model_ext.solve(
        analysis="kinematic",
        t_final=T_FINAL,
        t_eval=t_eval,
        ic_correct=IC_CORRECT,
    )

    # -- Mirror stroke: compute L_CYL_MIN from SimB --------------------------
    #
    #   When the right cyl extends to L_CYL the left follower is dragged to
    #   some L_CYL_MIN < L0_CYL (beyond its nominal stroke).  By the geometric
    #   mirror symmetry of the mechanism (both sides share identical local link
    #   lengths and angles; SHIFT only translates x), retracting the right cyl
    #   to that same L_CYL_MIN forces the left to extend to L_CYL — the exact
    #   kinematic mirror of Sim B.
    #
    #   L_CYL_MIN is computed from the yoke_left angle at t = T of Sim B using
    #   the crank-arm formula (ground pivot → rod_yoke tip distance).
    #
    _phi_L_end = uT_ext[-1, 3 * 0 + 2]           # yoke_left phi at t=T of Sim B
    _PIVOT_L   = np.array([-712.987, 142.348])    # yoke_left ground pivot [mm]
    _r_L       = np.hypot(36.005, 212.877)         # arm length: |pivot → rod_yoke| [mm]
    _a0_L      = np.arctan2(-212.877, 36.005)      # initial arm angle [rad]
    _rk_L      = _PIVOT_L + _r_L * np.array([
        np.cos(_a0_L + _phi_L_end),
        np.sin(_a0_L + _phi_L_end),
    ])                                              # rod_yoke global pos at end of Sim B
    L_CYL_MIN  = float(np.linalg.norm(_rk_L))     # barrel_left pivot is at (0, 0)
    DELTA_S_RET = L0_CYL - L_CYL_MIN
    print(f"\n[Full ROM]  Mirror stroke: L_CYL_MIN = {L_CYL_MIN:.3f} mm"
          f"  (Δ_ret = {DELTA_S_RET:.3f} mm  vs  Δ_ext = {DELTA_S:.3f} mm)")

    # -- Simulation A: rest -> maximum retraction (steer left) ---------------
    print("\n[Full ROM]  Sim A -- retraction  (L0_CYL -> L_CYL_MIN, steer left)")
    model_ret = build_model(f_start=L0_CYL, f_end=L_CYL_MIN)
    T_ret, uT_ret = model_ret.solve(
        analysis="kinematic",
        t_final=T_FINAL,
        t_eval=t_eval,
        ic_correct=IC_CORRECT,
    )

    # -- Concatenation: max_retraction -> rest -> max_extension --------------
    #
    #   Sim A reversed: reads from maximum retraction back to rest.
    #   The last step of Sim A reversed (= rest, t = 0) is dropped to avoid
    #   duplicating the boundary already present as step 0 of Sim B.
    #
    #   Combined time axis:  [-T_FINAL, ..., -dt, 0, dt, ..., T_FINAL]
    #

    # ── Diagnostic: check raw Sim A / Sim B final angles before concatenation
    _IL = 3 * 0 + 2   # yoke_left  phi column
    _IR = 3 * 4 + 2   # yoke_right phi column
    print("\n-- Diagnostic (raw solver output, before concatenation) " + "-" * 18)
    print("  Sim A (retraction, t=0→T):")
    print(f"    t=0 :  yoke_left={np.degrees(uT_ret[ 0, _IL]):+.3f}°  yoke_right={np.degrees(uT_ret[ 0, _IR]):+.3f}°  (should be ≈ 0)")
    print(f"    t=T :  yoke_left={np.degrees(uT_ret[-1, _IL]):+.3f}°  yoke_right={np.degrees(uT_ret[-1, _IR]):+.3f}°  (max retraction)")
    print("  Sim B (extension, t=0→T):")
    print(f"    t=0 :  yoke_left={np.degrees(uT_ext[ 0, _IL]):+.3f}°  yoke_right={np.degrees(uT_ext[ 0, _IR]):+.3f}°  (should be ≈ 0)")
    print(f"    t=T :  yoke_left={np.degrees(uT_ext[-1, _IL]):+.3f}°  yoke_right={np.degrees(uT_ext[-1, _IR]):+.3f}°  (max extension)")
    print("  Expected symmetry: SimA(t=T).left ≈ −SimB(t=T).right")
    print(f"    SimA left  = {np.degrees(uT_ret[-1, _IL]):+.3f}°  vs  −SimB right = {-np.degrees(uT_ext[-1, _IR]):+.3f}°")
    print(f"    SimA right = {np.degrees(uT_ret[-1, _IR]):+.3f}°  vs  −SimB left  = {-np.degrees(uT_ext[-1, _IL]):+.3f}°")
    print("-" * 72)
    # ── end diagnostic

    T_full  = np.concatenate([-T_ret[::-1][:-1],  T_ext])
    uT_full = np.concatenate([ uT_ret[::-1][:-1], uT_ext], axis=0)

    # Rebuild model_ext._result_container with the full 4001-step data so
    # that PostProcessor curve lengths match the T_full time vector.
    _patch_full_rom(model_ret, model_ext, T_full, uT_full)

    # -- ROM summary ---------------------------------------------------------
    #
    # Body order (from 'bodies' list):
    #   0 yoke_left | 1 barrel_left | 2 rod_left  | 3 wheel_left
    #   4 yoke_right| 5 barrel_right| 6 rod_right | 7 wheel_right | 8 tie_rod
    #
    # Orientation (body i):  uT[:, 3*i + 2]
    #
    idx_phi_yoke_left  = 3 * 0 + 2   #  2
    idx_phi_yoke_right = 3 * 4 + 2   # 14

    phi_yoke_left  = np.degrees(uT_full[:, idx_phi_yoke_left])
    phi_yoke_right = np.degrees(uT_full[:, idx_phi_yoke_right])

    # The rest configuration is at index N_EVAL - 1 (first point of Sim B).
    idx_rest = N_EVAL - 1
    phi_left_rest  = phi_yoke_left[idx_rest]
    phi_right_rest = phi_yoke_right[idx_rest]

    print()
    print("-- Full ROM summary " + "-" * 52)
    print(f"  Cylinder stroke:  [{-DELTA_S_RET:+.2f}, {+DELTA_S:+.2f}]  mm  (from rest, asymmetric)")
    print(f"  yoke_left   phi: [{phi_yoke_left.min():+.3f} deg, {phi_yoke_left.max():+.3f} deg]"
          f"  (rest = {phi_left_rest:+.3f} deg)")
    print(f"  yoke_right  phi: [{phi_yoke_right.min():+.3f} deg, {phi_yoke_right.max():+.3f} deg]"
          f"  (rest = {phi_right_rest:+.3f} deg)")
    print("-" * 72)

    # -- PostProcessor -------------------------------------------------------
    # model_ext is the reference model (same topology as model_ret).
    # Its _result_container has been patched with the full 4001-step data.
    # Time axis runs from -T_FINAL to +T_FINAL:
    #   negative half -> steering left  (retraction, Sim A reversed)
    #   positive half -> steering right (extension,  Sim B)
    PostProcessor(model=model_ext, T=T_full, uT=uT_full).show()
