# standard packages
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

#* developed functions
from functions import *

#! s_rot
a = np.array([2, 3])
b = np.array([[2], [3]])
validate_shape(b)
a_rotated = s_rot(a)
b_rotated = s_rot(b)

ecchime = 1 

#! A_matrix
theta = np.pi/6
s_local = np.array([[2], [3]])
A = A_matrix(theta)
s = A@s_local

#! example 3.1 
r = np.array([[2.5], [1.2]])
phi = np.radians(325)
s_A_local = np.array([[2.18], [0]])
s_B_local = np.array([[-1.8], [1.3]])
A = A_matrix(phi)
# a
s_A = A@s_A_local
s_B = A@s_B_local
# b
r_A = r_Point(r, s_A)
r_B = r_Point(r, s_B)
# c
s_BA = r_B - r_A

#! example 3.2
r_d = np.array([[1],[-2]])
phi_d = 1

r_dd = np.array([[-0.4],[0.4]])
phi_dd = 4.65 

r_A_d = r_Point_d(r_d, s_A, phi_d)
r_A_dd = r_Point_dd(r_dd, s_A, phi_d, phi_dd)

#! example 4.4
#* data
s_A_local, s_B_local = np.array([[0.15], [0]]), np.array([[0], [0.1]])
r1, r2 = np.array([[0.1], [0.05]]), np.array([[0.3], [-0.05]])
phi1, phi2 = np.pi/4, np.radians(20)
r1_d, r2_d = np.array([[0.1], [0.2]]), np.array([[-0.2], [0.1]])
phi1_d, phi2_d = -0.25, 0.12
L0, k, dc, fa = 0.15, 10, 5, -2

#* solution
A1  = A_matrix(phi1)
A2  = A_matrix(phi2)
r_A = my_r_Point(r1, s_A_local, A1)
r_B = my_r_Point(r2, s_B_local, A2)
d   = r_A - r_B

r_A_d = my_r_Point_d(r1_d, s_A_local, A1, phi1_d)
r_B_d = my_r_Point_d(r2_d, s_B_local, A2, phi2_d)
d_d   = r_A_d - r_B_d

f_sda = pp_sda(d, d_d, k, L0, dc, fa)
f_A = -f_sda 
f_B = f_sda

s_A = A1@s_A_local
s_B = A2@s_B_local
n_A = (s_rot(s_A).T)@f_A
n_B = (s_rot(s_B).T)@f_B

print(f_A)
print(f_B)
print(" ")
print(n_A)
print(n_B)

ecchime = 1