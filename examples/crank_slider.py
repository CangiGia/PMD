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
from pmd.core import *
from pmd.gui import PostProcessor


#%% Bodies
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


#%% Markers
# Ground markers
mk_g_crank = Ground.add_marker([0.0, 0.0])                       # crank pivot O1 = (0, 0)
mk_g_slider = Ground.add_marker([0.4644110316, 0.0], theta=0.0)  # slider guide (x-axis)

# Crankshaft markers (half-length = 0.05 m)
mk_c_ground = crankshaft.add_marker([-0.05, 0.0])  # ground-side end
mk_c_rod = crankshaft.add_marker([0.05, 0.0])      # rod-side end

# Rod markers (half-length = 0.2 m)
mk_r_crank = rod.add_marker([-0.2, 0.0])   # crank-side end
mk_r_slider = rod.add_marker([0.2, 0.0])   # slider-side end

# Slider markers (CM coincides with rod-slider joint)
mk_s_rod = slider.add_marker([0.0, 0.0])                # rev with rod
mk_s_ground = slider.add_marker([0.0, 0.0], theta=0.0)  # tran with ground


#%% Function (motion driver)
# phi_crank - phi_rod = c0 + c1*t
fn_mot = Function(
    type='a',
    coeff=[
        np.deg2rad(55.1821),  # c0: initial relative angle (phi_crank - phi_rod) [rad]
        np.deg2rad(20.0),     # c1: angular velocity [rad/s]
        0.0,                  # c2
    ],
)


#%% Joints
j_gc = RevJoint(
    iMarker=mk_g_crank,
    jMarker=mk_c_ground,
    name='ground_crankshaft_joint',
)

j_cr = RevJoint(
    iMarker=mk_c_rod,
    jMarker=mk_r_crank,
    name='crankshaft_rod_joint',
)

j_rs = RevJoint(
    iMarker=mk_s_rod,
    jMarker=mk_r_slider,
    name='rod_slider_joint',
)

j_gs = TranJoint(
    iMarker=mk_s_ground,
    jMarker=mk_g_slider,
    name='ground_slider_joint',
)

j_mot = RelRotJoint(
    iBody=crankshaft,
    jBody=rod,
    iFunct=fn_mot,
    name='crank_rod_motor',
)


#%% Forces
fw = Weight(gravity=9.80665)  # m/s^2, along -y


#%% Simulation
model = PlanarMultibodyModel(
    bodies=[crankshaft, rod, slider],
    joints=[j_gc, j_cr, j_rs, j_gs, j_mot],
    forces=[fw],
    functions=[fn_mot],
)

T, uT = model.solve(
    method='CASADI-DAE',
    t_final=20,
    t_eval=np.linspace(0, 20, 20000),
    ic_correct=True,
)

post_proc = PostProcessor(model=model, T=T, uT=uT)
post_proc.show()