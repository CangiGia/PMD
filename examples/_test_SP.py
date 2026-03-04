import numpy as np
import scipy as sc
from PMD.src.builder import *
from PMD.src.mechanics import *
from PMD.src.solver import *
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

#%% points ...
p0 = Point() #// 0 is always the ground ...
p0.body = Ground
p0.sPlocal = np.array([0., 0.2])

p1 = Point()
p1.body = b1
p1.sPlocal = np.array([0., 0.])

p2 = Point()
p2.body = b2
p2.sPlocal = np.array([0., 0.5])

#%% unit vectors ...
u0 = uVector()
u0.body = Ground
u0.ulocal = np.array([1., 0.]) # parallel to the x-axis

u1 = uVector()
u1.body = b1
u1.ulocal = np.array([1., 0.]) # parallel to the x-axis

#%% forces ...
f0_1 = Force()
f0_1.type = "ptp"
f0_1.iPoint = p1
f0_1.jPoint = p0
f0_1.k = 20.0
f0_1.L0 = 0.6

fw = Force()
fw.type = "weight"

#%% joints ...
j0_1 = Joint()
j0_1.type = "tran"
j0_1.iPoint = p1
j0_1.jPoint = p0
j0_1.iUvec = u1
j0_1.jUvec = u0

j1_2 = Joint()
j1_2.type = "rev"
j1_2.iPoint = p1
j1_2.jPoint = p2

#%% model simulation
my_dynamic_model = PlanarMultibodyModel(
    bodies=[b1, b2],
    joints=[j0_1, j1_2],
    forces=[f0_1, fw],
    points=[p0, p1, p2],
    uvectors=[u0, u1])
T, uT = my_dynamic_model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='SP.txt', model_title='SP')
