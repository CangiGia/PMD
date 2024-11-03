import numpy as np
from scipy.optimize import fsolve
from PMD_functions import *

#* Data
L = 0.5
m1, J1 = 2, 0.04
m2, J2 = 1, 0.1
kr = 6
theta0, theta1 = 0.0, 60
k, L0, dc = 20, 0.7, 5
dAO_local = 0.8
g = 9.81

#* Solution
A = A_matrix(np.radians(theta1))
phi1, dphi1 = np.pi/3, 0
phi2, dphi2 = phi1, 0 
sO1_local = colvect(0, 0.5)
sA2_local = colvect(0, 0)
r1, dr1 = colvect(0.5*np.sin(np.radians(theta1)), -0.5*np.cos(np.radians(theta1))), colvect(0,0)
r2, dr2 = colvect(0.8*np.sin(np.radians(theta1)), -0.8*np.cos(np.radians(theta1))), colvect(0,0)
uxi, ueta = colvect(1, 0), colvect(0, 1)

# Mass matrix
mdiag = [m1, m1, J1, m2, m2, J2]
M = np.diag(mdiag)

# Jacobian matrix for a revolute joint
sO1 = A@sO1_local
Dr = np.hstack([np.eye(2), s_rot(sO1)])

# Jacobian matrix for a traslational joint
u1, u2 = A@ueta, A@uxi
d = r2 - r1
Dt1 = np.hstack([np.vstack([-s_rot(u1).T, np.zeros([1,2])]), np.vstack([-u1.T@d, -1])])
Dt2 = np.hstack([np.vstack([s_rot(u1).T, np.zeros([1,2])]), np.vstack([0, 1])])

# complete Jacobian matrix
D = np.hstack([np.vstack([Dr, Dt1]), np.vstack([np.zeros([2,3]), Dt2])])

# weight forces
uy = colvect(0, 1)
fw1, fw2 = -m1*g*uy, -m2*g*uy

# spring-damper element forces
fsd = pp_sd(r2, dr2, k, L0, dc)
fsd2 = -fsd

# tosrional spring moment
mts = r_s(np.radians(theta1), kr, theta0)
t1 = -mts

# matrix coefficient of constraint accelerations
gammar = (sO1*dphi1**2)
gammat = np.vstack([((s_rot(u1).T@d)*dphi1 + 2*(u1.T@dr1))*dphi1, np.zeros([1,1])])
gamma = np.vstack([gammar, gammat])

# EQMs definition 
lhs = np.hstack([np.vstack([M, D]), np.vstack([D.T, np.zeros([4,4])])])
rhs = np.vstack([fw1, t1, fw2+fsd2, 0, gamma])

solutions = np.linalg.solve(lhs, rhs)
dcc = solutions[0:6]
lm = solutions[6:10]

ecchime = 1 