import numpy as np
import scipy as sc
from maker import *
from solver import *
from functions import *
import matplotlib.pyplot as plt


#* Multi-Body model creation - Bodies, Joints, Forces
#%% bodies
b1 = Body()
b1.m = 2
b1.J = 0.5
b1.r = np.array([0.4398, 0.2512])
b1.p = -0.0367

b2 = Body()
b2.m = 30
b2.J = 2.5 
b2.r = np.array([0.6817, 0.3498])
b2.p = 0.0783

b3 = Body()
b3.m = 1
b3.J = 0.5 
b3.r = np.array([0.4463, 0.4308])
b3.p = 6.5222

#%% points
f0 = Point()
f0.Bindex = 0
f0.sPlocal = np.array([0.38, 0.43])

o0 = Point()
o0.Bindex = 0
o0.sPlocal = np.array([0.32, 0.40])

q0 = Point()
q0.Bindex = 0
q0.sPlocal = np.array([0.20, 0.26])

q1 = Point()
q1.Bindex = 1
q1.sPlocal = np.array([-0.24, 0.0])

a1 = Point()
a1.Bindex = 1
a1.sPlocal = np.array([0.18, 0.0])

e1 = Point()
e1.Bindex = 1
e1.sPlocal = np.array([0.0, 0.0])

a2 = Point()
a2.Bindex = 2
a2.sPlocal = np.array([-0.07, -0.10])

b2 = Point()
b2.Bindex = 2
b2.sPlocal = np.array([-0.10, 0.12])

b3 = Point()
b3.Bindex = 3
b3.sPlocal = np.array([0.13, 0.0])

o3 = Point()
o3.Bindex = 3
o3.sPlocal = np.array([-0.13, 0.0])

#%% joints
j1 = Joint() #! Revolute joint in Q 
j1.type = 'rev'
j1.iPindex = 1
j1.jPindex = 8

j2 = Joint() #! Revolute joint in A 
j2.type = 'rev'
j2.iPindex = 2
j2.jPindex = 3

j3 = Joint() #! Revolute joint in B 
j3.type = 'rev'
j3.iPindex = 4
j3.jPindex = 5

j4 = Joint() #! Revolute joint in O 
j4.type = 'rev'
j4.iPindex = 6
j4.jPindex = 7

ecchime = 1