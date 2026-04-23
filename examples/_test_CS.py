"""
CS (Crank-Slider) Model - Planar Multibody Dynamics
=====================================================
Translated from Adams/View model: PMD validation/crank_slider*.cmd

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
from PMD.src import *
from gui.app import PostProcessor

#%% Bodies
crankshaft = Body(
    mass=0.8169313258, # kg
    inertia=1.294122161962e-4, # kg*m^2
    position=[0.03535533906, 0.03535533906], # m
    orientation=np.deg2rad(45.0),
    name='crankshaft',
)

rod = Body(
    mass=2.689171326, # kg
    inertia=4.414522158e-4, # kg*m^2
    position=[0.2675608549, 0.03535522707], # m
    orientation=np.deg2rad(-10.1821),
    name='rod',
)

slider = Body(
    mass=0.2614140191, # kg
    inertia=4.182624305e-5, # kg*m^2
    position=[0.4644110316, 0.0], # m
    name='slider',
) 

#%% Markers
# Ground markers
mk_g_crank  = Ground.add_marker([0.0, 0.0])                        # rev with crankshaft
mk_g_slider = Ground.add_marker([0.4644110316, 0.0], theta=0.0)   # tran with slider (x-dir)

# Crankshaft markers  (half-length = 0.05 m)
mk_c_ground = crankshaft.add_marker([-0.05, 0.0])  # at global (0, 0)
mk_c_rod    = crankshaft.add_marker([ 0.05, 0.0])  # at crankshaft tip

# Rod markers  (half-length = 0.2 m)
mk_r_crank  = rod.add_marker([-0.2, 0.0])          # at crank-rod joint
mk_r_slider = rod.add_marker([ 0.2, 0.0])          # at rod-slider joint

# Slider markers  (CM coincides with rod-slider joint)
mk_s_rod    = slider.add_marker([0.0, 0.0])                       # rev with rod
mk_s_ground = slider.add_marker([0.0, 0.0], theta=0.0)            # tran with ground

#%% Function  (motion driver: 20 deg/s constant angular velocity)
# RelRotJoint(iBody=crankshaft, jBody=rod):  phi_crank - phi_rod = f(t)
#   f(t) = c0 + c1*t  with c0 = initial relative angle = 55.1821 deg in rad
fn_mot = Function(
    type='a',
    coeff=[np.deg2rad(55.1821),   # c0: initial relative angle (phi_crank - phi_rod)
           np.deg2rad(20.0),      # c1: 20 deg/s
           0.0],                  # c2
)

#%% Joints
j_gc  = RevJoint(iMarker=mk_g_crank,  jMarker=mk_c_ground, name='ground_crankshaft_joint')
j_cr  = RevJoint(iMarker=mk_c_rod,    jMarker=mk_r_crank,  name='crankshaft_rod_joint')
j_rs  = RevJoint(iMarker=mk_s_rod,    jMarker=mk_r_slider, name='rod_slider_joint')
j_gs  = TranJoint(iMarker=mk_s_ground,jMarker=mk_g_slider, name='ground_slider_joint')
j_mot = RelRotJoint(iBody=crankshaft, jBody=rod, iFunct=fn_mot, name='crank_rod_motor')

#%% Forces
fw = Weight(gravity=9.80665)

#%% Simulation

model = PlanarMultibodyModel(
    bodies=[crankshaft, rod, slider],
    joints=[j_gc, j_cr, j_rs, j_gs, j_mot],
    forces=[fw],
    functions=[fn_mot],
)

T, uT = model.solve(
    method='Radau',
    t_final=5,
    t_eval=np.linspace(0, 5, 5000),
    ic_correct=True,
    baumgarte_alpha=0.0,
    baumgarte_beta=0.0,
)

# if __name__ == '__main__':
#     from PMD.examples._plot_utils import plot_comparison
#     plot_comparison(T, uT, matlab_filename='CS.txt', model_title='CS')

post_proc = PostProcessor(model=model, T=T, uT=uT)
post_proc.show()