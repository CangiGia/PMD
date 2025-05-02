import numpy as np
import scipy as sc
from PMD.src.builder import *
from PMD.src.solver import *
from PMD.src.functions import *
import matplotlib.pyplot as plt


#* Multi-Body model creation - Bodies, Joints, Forces
#%% bodies
B1 = Body() # lower suspension arm
B1.m = 2
B1.J = 0.5
B1.r = np.array([0.4398, 0.2512])
B1.p = 6.2465

B2 = Body() # wheel assembly
B2.m = 30
B2.J = 2.5 
B2.r = np.array([0.6817, 0.3498])
B2.p = 0.0783

B3 = Body() # upper suspension arm
B3.m = 1
B3.J = 0.5 
B3.r = np.array([0.4463, 0.4308])
B3.p = 6.5222

#%% points
q1 = Point() # Q point - lower suspension arm side, body 1
q1.Bindex = 1
q1.sPlocal = np.array([-0.24, 0.0])

a1 = Point() # A point - lower suspension arm side, body 1
a1.Bindex = 1
a1.sPlocal = np.array([0.18, 0.0])

a2 = Point() # A point - wheel assembly side, body 2
a2.Bindex = 2
a2.sPlocal = np.array([-0.07, -0.10])

b2 = Point() # B point - wheel assembly side, body 2
b2.Bindex = 2
b2.sPlocal = np.array([-0.10, 0.12])

b3 = Point() # B point - upper suspension arm side, body 3
b3.Bindex = 3
b3.sPlocal = np.array([0.13, 0.0])

o3 = Point() # O point - upper suspension arm side, body 3
o3.Bindex = 3
o3.sPlocal = np.array([-0.13, 0.0])

o0 = Point() # O point - ground
o0.Bindex = 0
o0.sPlocal = np.array([0.32, 0.40])

q0 = Point() # Q point - ground
q0.Bindex = 0
q0.sPlocal = np.array([0.20, 0.26])

e1 = Point() # E point - Lower suspension arm side, body 1
e1.Bindex = 1 
e1.sPlocal = np.array([0.0, 0.0])

f0 = Point() # F point - Ground
f0.Bindex = 0
f0.sPlocal = np.array([0.38, 0.43])

#%% joints
j1 = Joint() #! Revolute joint in Q 
j1.type = 'rev'
j1.iPindex = 0
j1.jPindex = 7

j2 = Joint() #! Revolute joint in A 
j2.type = 'rev'
j2.iPindex = 1
j2.jPindex = 2

j3 = Joint() #! Revolute joint in B 
j3.type = 'rev'
j3.iPindex = 3
j3.jPindex = 4

j4 = Joint() #! Revolute joint in O 
j4.type = 'rev'
j4.iPindex = 5
j4.jPindex = 6

#%% forces - custom forces
s1 = Force()
s1.type = 'ptp'
s1.iPindex = 8
s1.jPindex = 9
s1.k = 90000
s1.L0 = 0.23
s1.dc = 1100

# custom mathematical model
def my_force(B2, s2):
    """
    Custom force used to define the wheel contact condition.
    """

    # unilateral spring-damper in the y-direction
    dely = B2.r[1] - s2.L0
    if dely < 0:
        fy = (s2.k * dely + s2.dc * B2.dr[1]).item()
        fsd = np.array([0, -fy])
        B2._f = B2._f + fsd.reshape(2,1)

s2 = Force() # wheel contact force
s2.type = 'user'
s2.k = 50000
s2.L0 = 0.35
s2.dc = 1000
s2.callback = my_force

s3 = Force()
s3.type = 'weight'
s3.wgt = np.array([0, -1])

#%% solution
my_dynamic_model = PlanarMultibodyModel(verbose=True)
time, solution = my_dynamic_model.solve(method='RK45')

plt.figure()
plt.plot(time, solution[:,4])
plt.show()

# np.savetxt('uT_python.txt', solution)
# np.savetxt('T_python.txt', time)

ecchime = 1