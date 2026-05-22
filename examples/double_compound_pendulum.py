""" DCP (Double Compound Pendulum) Model - Planar Multibody Dynamics (pmd)

      (O1)
        *
         * 
          *         
           * 
           (O2) * * * * * * (O3)
                 
                
BODY PARAMETERS
---------------
+--------+-----------+------------------+---------------------+------------------+
| Body   | Mass [kg] | Inertia [kg·m²]  | IC: CM pos [m]      | IC: angle [deg]  |
+--------+-----------+------------------+---------------------+------------------+
| link_1 | 1.756261  |    0.011854      | (0.1218, −0.0281)   | −12.9946 (=347°) |
| link_2 | 2.692381  |    0.042171      | (0.4378, −0.0086)   | +13.7698         |
+--------+-----------+------------------+---------------------+------------------+
"""

import numpy as np
from pmd.core import (
    Body,
    Ground,
    PlanarMultibodyModel,
    RevJoint,
    Weight,
)
from pmd.core.shapes import Link


T_FINAL    = 35.0
N_EVAL     = 35001
IC_CORRECT = True

_L1 = 250.0e-3   # m
_L2 = 400.0e-3   # m
_W  =  40.0e-3   # m
_D  =  20.0e-3   # m

_THETA1_DEG = -12.9946167919 
_THETA2_DEG =  13.7697807428 

_M1 = 1.756261   # kg
_J1 = 0.011854   # kg·m²
_M2 = 2.692381   # kg
_J2 = 0.042171   # kg·m²

_th1 = np.deg2rad(_THETA1_DEG)
_th2 = np.deg2rad(_THETA2_DEG)

link_1 = Body(
    mass=_M1,
    inertia=_J1,
    position=[_L1 / 2 * np.cos(_th1), _L1 / 2 * np.sin(_th1)],
    orientation=_th1,
    name='link_1',
    shape=Link(length=_L1, thickness=_W),
)

link_2 = Body(
    mass=_M2,
    inertia=_J2,
    position=[_L1 * np.cos(_th1) + _L2 / 2 * np.cos(_th2),
              _L1 * np.sin(_th1) + _L2 / 2 * np.sin(_th2)],
    orientation=_th2,
    name='link_2',
    shape=Link(length=_L2, thickness=_W),
)

g_l1  = Ground.add_marker([0.0, 0.0])         # O1 - ground side 
l1_g  = link_1.add_marker([-_L1 / 2, 0.0])    # O1 - link_1 side
l1_l2 = link_1.add_marker([ _L1 / 2, 0.0])    # O2 - link_1 side
l2_l1 = link_2.add_marker([-_L2 / 2, 0.0])    # O2 - link_2 side

j_ground_l1 = RevJoint(iMarker=g_l1,  jMarker=l1_g,  q0=_th1,        name='ground_link_1_joint')
j_l1_l2     = RevJoint(iMarker=l1_l2, jMarker=l2_l1, q0=_th2 - _th1, name='link_1_link_2_joint')

fw = Weight(gravity=9.80665)

model = PlanarMultibodyModel(
    bodies=[link_1, link_2],
    joints=[j_ground_l1, j_l1_l2],
    forces=[fw],
)


if __name__ == "__main__":
    from pmd.gui import PostProcessor

    T, uT = model.solve(
        t_final=T_FINAL,
        t_eval=np.linspace(0, T_FINAL, N_EVAL),
        ic_correct=IC_CORRECT,
    )
    PostProcessor(model=model, T=T, uT=uT).show()
