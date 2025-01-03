#! This example is the continue of the Example 6_1
import numpy as np
from scipy.optimize import fsolve
from functions import * 


#* Data
m1, J1 = 2, 1.5
m2, J2 = 1.5, 0.5
k, L0 = 100, 2
dc = 20
g = 9.81 
r1A, r1B, r1G = colvect(-1.5, 3), colvect(-1.5, 0), colvect(-1.5, 1.5)
r2C, r2G = colvect(-1.5, 4.5), colvect(-0.5, 4.5)
dr1B = colvect(-2, 1)
ux, uy = colvect(1, 0), colvect(0, 1)
w1, w2 = -m1*g*uy, -m2*g*uy

#* Solution 
# spring force 
ds = r2C - r1A
fs = pp_s(ds, k, L0) # function to retrive the spring force
fsA, fsC = fs*uy, -fs*uy

# damper force
dd = r1B
fd = pp_sd(dd, dr1B, 0, 0, dc) # function to retrieve the damper force
fdB = -fd*ux

# moments evaluation and generalized force vectors definition 
s1A = r1A - r1G
s1B = r1B - r1G
s2C = r2C - r2G
h1 = np.vstack((w1+fsA+fdB, s_rot(s1A).T@fsA + s_rot(s1B).T@fdB))
h2 = np.vstack((w2+fsC, s_rot(s2C).T@fsC))

# mass matrix evaluation
diag = [m1, m1, J1, m2, m2, J2]
M = np.diag(diag)

# +++++++++++
# Contiue 
# +++++++++++
r2D = colvect(0.5, 4.5)
r2D_hat = unit_vector(r2D) # force components due to the reaction force
s2D = r2D - r2G 
n2 = s_rot(s2D).T@r2D_hat # moment componet due to the reaction force

ecchime = 1