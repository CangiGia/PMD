import numpy as np
import scipy as sc
from src.maker import *
from src.solver import *
from src.functions import *
import matplotlib.pyplot as plt


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
p0.Bindex = 0
p0.sPlocal = np.array([0., 0.2])

p1 = Point()
p1.Bindex = 1
p1.sPlocal = np.array([0., 0.])

p2 = Point()
p2.Bindex = 2
p2.sPlocal = np.array([0., 0.5])

#%% unit vectors ...
u0 = uVector()
u0.Bindex = 0 
u0.ulocal = np.array([1., 0.]) # parallel to the x-axis

u1 = uVector()
u1.Bindex = 1 
u1.ulocal = np.array([1., 0.]) # parallel to the x-axis

#%% forces ...
f0_1 = Force()
f0_1.type = "ptp"
f0_1.iPindex = 1
f0_1.jPindex = 0
f0_1.k = 20.0
f0_1.L0 = 0.6

fw = Force()
fw.type = "weight"

#%% joints ...
j0_1 = Joint()
j0_1.type = "tran"
j0_1.iPindex = 1
j0_1.jPindex = 0
j0_1.iUindex = 1
j0_1.jUindex = 0

j1_2 = Joint()
j1_2.type = "rev"
j1_2.iPindex = 1
j1_2.jPindex = 2

#%% some controls ... 
body_count    = Body.get_count()
point_count   = Point.get_count()
uvector_count = uVector.get_count()
force_count   = Force.get_count()
joint_count   = Joint.get_count()

print(f" ")
print(f"\t ... some controls ...")
print(f"\t ... number of Body instances:    {body_count}   ")
print(f"\t ... number of Point instances:   {point_count}  ")
print(f"\t ... number of uVector instances: {uvector_count}")
print(f"\t ... number of Force instances:   {force_count}  ")
print(f"\t ... number of Joint instances:   {joint_count}  ")

#%% model simulation ...
my_dynamic_model = PlanarDynamicModel(verbose=True)
time, solution = my_dynamic_model.solve(method='RK45')

plt.figure()
plt.plot(time, solution[:,4])
plt.show()

np.savetxt('uT_python.txt', solution)
np.savetxt('T_python.txt', time)

ecchime = 1