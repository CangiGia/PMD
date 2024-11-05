# continue of example 7.1
import numpy as np
from scipy.optimize import fsolve
from PMD.pmd_functions import * 


#* Data
g = 9.81
l1, l2 = 0.3, 0.24
m1, J1 = 0.2, 0.03
m2, J2 = 0.15, 0.02
k, l0, dc = 50, 0.2, 20
r1, dr1 = colvect(0.02, 0.2), colvect(-0.05, 0.1)
r2, dr2 = colvect(0.22, 0.1), colvect(-0.09, 0.13)
phi1, dphi1 = np.radians(45), -0.3
phi2, dphi2 = np.radians(15), 0.07

#* Solution
mdiag = (m1,m1,J1,m2,m2,J2)
M = np.diag(mdiag)
sA_local = colvect(l1/2, 0)
sB_local = colvect(-l1/2, 0)
sC_local = colvect(0, l2/2)
sD_local = colvect(0, -l2/2)

# weight force
ux = colvect(1, 0)
uy = colvect(0, 1)
fw1, fw2 = -m1*g*uy, -m2*g*uy

# elastic force
A1, A2 = A_matrix(phi1), A_matrix(phi2)
rA = my_r_Point(r1, sA_local, A1)
rC = my_r_Point(r2, sC_local, A2)
ds = rA - rC
fs = pp_s(ds, k, l0)

# damping force
rB = my_r_Point(r1, sB_local, A1)
rQ = colvect(0, 0)
dd = rB - rQ 
drQ = colvect(0, 0)
drB = my_r_Point_d(dr1, sB_local, A1, dphi1)
ddd = drB - drQ
fd = pp_sd(dd, ddd, k = 0, L0 = l0, dc = dc)

# moments due to active forces
sA = A1@sA_local
sB = A1@sB_local
sC = A2@sC_local
fAs = -fs
fCs = fs
fBd = -fd
n1 = (s_rot(sA).T@fAs) + (s_rot(sB).T@fBd)
n2 = (s_rot(sC).T@fCs)

# rhsa of the system
rhsa = np.vstack([(fAs+fBd+fw1), n1, (fCs+fw2), n2])

# computing acceleration
ddc = np.linalg.solve(M,rhsa)
print("Computed accelerations array: \n")
print(ddc)

##### ##### ##### ##### #####
# Adding a rigid rod between point B and D
# Determining the Lagrangian multiplier
# and the Jacobian matrix of the system
##### ##### ##### ##### #####
Lrod = 0.355
sD = A2@sD_local
rD = my_r_Point(r2, sD_local, A2)
drD = my_r_Point_d(rD, sD_local, A2, dphi2)
drod = rB - rD
urod = drod/Lrod
ddrod = drB - drD
durod = ddrod/Lrod

# check on posistioin constraint evaluation
PHI = 0.5*(urod.T@drod - Lrod)
print(' \n')
print('Phi value: ', PHI)
print(' \n')

# Jacobian matrix definition
D_rev_rev = np.hstack([urod.T, urod.T@s_rot(sB), -urod.T, -urod.T@s_rot(sD)])

# velocity vector evaluation, check on velocity violation
cd = np.vstack([np.vstack([dr1, dphi1]), np.vstack([dr2, dphi2])])
dPhi = D_rev_rev@cd
print('Phi_d value: ', dPhi)
print(' \n')

# gamma evaluation to define the rhs of the aceelerations constraint
dsB = (s_rot(sB)*dphi1)
dsD = (s_rot(sD)*dphi2)
gamma = (-durod.T@ddrod)+(s_rot(urod).T@((dsB*dphi1)-(dsD*phi2)))
print('gamma: ', gamma)
print(' \n')

# solution 
lhs = np.hstack([np.vstack([M, D_rev_rev]), np.vstack([-D_rev_rev.T,0])])
rhs = np.vstack([rhsa, gamma])
solution = np.linalg.solve(lhs, rhs)
print("------------\n")
print("Solution:\n")
print("------------\n")
print(solution)

ecchime = 1