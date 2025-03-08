import numpy as np
import scipy as sc
from PMD.maker import *
from PMD.solver import *
from PMD.functions import *
import matplotlib.pyplot as plt


#* Multi-Body model creation - Bodies, Joints, Forces
# bodies
b1 = Body()
b1.m = 5
b1.J = 4
b1.r = np.array([1, 0.2])

b2 = Body()
b2.m = 2
b2.J = 0.2 
b2.r = np.array([1.25, -0.233])
b2.p = np.pi/6

# points
p0 = Point()
p0.Bindex = 0
p0.sPlocal = np.array([0, 0.2])

p1 = Point()
p1.Bindex = 1
p1.sPlocal = np.array([0, 0])

p2 = Point()
p2.Bindex = 2
p2.sPlocal = np.array([0, 0.5])

# unit vectors
u1 = uVector()
u1.Bindex = 0 
u1.ulocal = np.array([1, 0])

u2 = uVector()
u2.Bindex = 1
u2.ulocal = np.array([1, 0])

# forces
f1 = Force()
f1.type = "ptp"
f1.iPindex = 1
f1.jPindex = 0
f1.k = 20
f1.L0 = 0.6

f2 = Force()
f2.type = "weight"
f2.gravity = 9.81
f2.wgt = np.array([0, -1])

# joints
j1 = Joint()
j1.type = "tran"
j1.iPindex = 1
j1.jPindex = 0
j1.iUindex = 1
j1.jUindex = 0

j2 = Joint()
j2.type = "rev"
j2.iPindex = 1
j2.jPindex = 2

# // ... some controls ... 
body_count    = Body.get_count()
point_count   = Point.get_count()
uvector_count = uVector.get_count()
force_count   = Force.get_count()
joint_count   = Joint.get_count()

print(f" ")
print(f"\t ... some controls ...")
print(f"\t ... number of Body instances:    {body_count}    ...")
print(f"\t ... number of Point instances:   {point_count}   ...")
print(f"\t ... number of uVector instances: {uvector_count} ...")
print(f"\t ... number of Force instances:   {force_count}   ...")
print(f"\t ... number of Joint instances:   {joint_count}   ...")

my_dynamic_model = PlanarDynamicModel(_verbose=True)
time, solution = my_dynamic_model.solve()

plt.figure()
plt.plot(time, solution[:,5])
plt.show() 

np.savetxt('uT_python.txt', solution)
np.savetxt('T_python.txt', time)

ecchime = 1