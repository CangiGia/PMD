"""
CS (Crank-Slider) Model - Planar Multibody Dynamics
=====================================================
Translated from Adams/View model: pmd validation/crank_slider*.cmd

Crank-slider mechanism: a crankshaft pinned to the ground drives a
connecting rod attached to a slider that translates along the x-axis.
A rotational motion driver on the crankshaft-rod joint prescribes the
relative rotation at 20 deg/s (constant angular velocity).

Bodies:
  crankshaft: link 0.1 m long, 45 deg from horizontal
  rod:        link 0.4 m long, -10.1821 deg from horizontal
  slider:     sphere R = 0.02 m, sliding along the x-axis

Joints:
  j_gc : revolute       Ground     <-> Crankshaft  at (0, 0)
  j_cr : revolute       Crankshaft <-> Rod         at crankshaft tip
  j_rs : revolute       Rod        <-> Slider      at rod tip
  j_gs : translational  Ground     <-> Slider      along x-axis
  j_mot: rel-rot        Crankshaft relative to Rod, driven by fn_mot

Forces:
  fw: weight (gravity along -y)
"""

import numpy as np
from pmd.core import *
from pmd.core.model import _GroundType
from pmd.gui import PostProcessor

# ---------------------------------------------------------------------------
# Factory — rebuilds the full CS model from scratch each call.
# Resets the Ground singleton markers and the Body counter so that
# successive calls produce independent, non-aliased instances.
# ---------------------------------------------------------------------------

def _make_cs_model():
    """Return a fresh PlanarMultibodyModel for the crank-slider mechanism."""
    _GroundType._instance._markers = [_GroundType._instance.origin]
    Body.COUNT = 0

    crankshaft = Body(
        mass=0.8169313258,
        inertia=1.294122161962e-4,
        position=[0.03535533906, 0.03535533906],
        orientation=np.deg2rad(45.0),
        name='crankshaft',
    )
    rod = Body(
        mass=2.689171326,
        inertia=4.414522158e-4,
        position=[0.2675608549, 0.03535522707],
        orientation=np.deg2rad(-10.1821),
        name='rod',
    )
    slider = Body(
        mass=0.2614140191,
        inertia=4.182624305e-5,
        position=[0.4644110316, 0.0],
        name='slider',
    )

    # Ground markers
    mk_g_crank  = Ground.add_marker([0.0, 0.0])
    mk_g_slider = Ground.add_marker([0.4644110316, 0.0], theta=0.0)

    # Crankshaft markers (half-length = 0.05 m)
    mk_c_ground = crankshaft.add_marker([-0.05, 0.0])
    mk_c_rod    = crankshaft.add_marker([ 0.05, 0.0])

    # Rod markers (half-length = 0.2 m)
    mk_r_crank  = rod.add_marker([-0.2, 0.0])
    mk_r_slider = rod.add_marker([ 0.2, 0.0])

    # Slider markers
    mk_s_rod    = slider.add_marker([0.0, 0.0])
    mk_s_ground = slider.add_marker([0.0, 0.0], theta=0.0)

    # Motion driver: phi_crank - phi_rod = c0 + c1*t  (20 deg/s)
    fn_mot = Function(
        type='a',
        coeff=[np.deg2rad(55.1821), np.deg2rad(20.0), 0.0],
    )

    j_gc  = RevJoint(iMarker=mk_g_crank,  jMarker=mk_c_ground, name='ground_crankshaft_joint')
    j_cr  = RevJoint(iMarker=mk_c_rod,    jMarker=mk_r_crank,  name='crankshaft_rod_joint')
    j_rs  = RevJoint(iMarker=mk_s_rod,    jMarker=mk_r_slider, name='rod_slider_joint')
    j_gs  = TranJoint(iMarker=mk_s_ground,jMarker=mk_g_slider, name='ground_slider_joint')
    j_mot = RelRotJoint(iBody=crankshaft, jBody=rod, iFunct=fn_mot, name='crank_rod_motor')

    fw = Weight(gravity=9.80665)

    return PlanarMultibodyModel(
        bodies=[crankshaft, rod, slider],
        joints=[j_gc, j_cr, j_rs, j_gs, j_mot],
        forces=[fw],
        functions=[fn_mot],
    )


# ---------------------------------------------------------------------------
# Simulation runs — each on an independent model instance
# ---------------------------------------------------------------------------

T_EVAL = np.linspace(0, 5, 5000)

#%% Run 1 — CasADi DAE (Radau collocation, constraints enforced exactly)
model_dae = _make_cs_model()
T_DAE, uT_DAE = model_dae.solve(
    method='CASADI-DAE',
    t_eval=T_EVAL,
    ic_correct=True,
)

#%% Run 2 — ODE / Radau, no stabilisation (pure constraint drift)
model_radau = _make_cs_model()
T_Radau, uT_Radau = model_radau.solve(
    method='Radau',
    t_eval=T_EVAL,
    ic_correct=True,
    baumgarte_alpha=0.0,
    baumgarte_beta=0.0,
)

#%% Run 3 — ODE / Radau + Baumgarte stabilisation
model_baum = _make_cs_model()
T_Baum, uT_Baum = model_baum.solve(
    method='Radau',
    t_eval=T_EVAL,
    ic_correct=True,
    baumgarte_alpha=5.0,
    baumgarte_beta=5.0,
)

# ---------------------------------------------------------------------------
# PostProcessor — all three sessions loaded; switch between them via
# Simulations menu in the left sidebar.
# ---------------------------------------------------------------------------

post_proc = PostProcessor(model=model_dae,   T=T_DAE,   uT=uT_DAE,   name="CASADI-DAE")
post_proc.add(            model_radau,        T_Radau,   uT_Radau,    name="Radau (no stabilisation)")
post_proc.add(            model_baum,         T_Baum,    uT_Baum,     name="Radau + Baumgarte (α=β=5)")
post_proc.show()

# ---------------------------------------------------------------------------
# Optional: ADAMS comparison (uncomment to run)
# ---------------------------------------------------------------------------
# from examples._plot_utils import plot_comparison
# plot_comparison(T_DAE, uT_DAE, matlab_filename='CS.txt', model_title='CS')