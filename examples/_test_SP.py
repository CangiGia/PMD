import numpy as np
import scipy as sc
from PMD.src import *
import matplotlib.pyplot as plt
from PMD.examples._plot_utils import plot_comparison


#* Multi-Body model creation - Bodies, Joints, Forces
#%% bodies ...
b1 = Body()
b1.m = 5.
b1.J = 4.
b1.r = np.array([1., 0.2])

b2 = Body()
b2.m = 2.
b2.J = 0.2
b2.r = np.array([1.25, -0.233])
b2.p = np.pi/6

#%% markers (merge Point + uVector into single Marker with theta for tran joint)
# p0 was Point(Ground, [0,0.2]) and u0 was uVector(Ground, [1,0]) -> theta=0.0
p0 = Ground.add_marker([0., 0.2], theta=0.0)

# p1 was Point(b1, [0,0]) and u1 was uVector(b1, [1,0]) -> theta=0.0
p1 = b1.add_marker([0., 0.], theta=0.0)

p2 = b2.add_marker([0., 0.5])

#%% forces ...
f0_1 = Force()
f0_1.type = "ptp"
f0_1.iMarker = p1
f0_1.jMarker = p0
f0_1.k = 20.0
f0_1.L0 = 0.6

fw = Force()
fw.type = "weight"

#%% joints ...
j0_1 = Joint()
j0_1.type = "tran"
j0_1.iMarker = p1
j0_1.jMarker = p0

j1_2 = Joint()
j1_2.type = "rev"
j1_2.iMarker = p1
j1_2.jMarker = p2

#%% model simulation
my_dynamic_model = PlanarMultibodyModel(
    bodies=[b1, b2],
    joints=[j0_1, j1_2],
    forces=[f0_1, fw])
T, uT = my_dynamic_model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='SP.txt', model_title='SP')
