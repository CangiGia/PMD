import numpy as np
from scipy.optimize import fsolve
from PMD.functions import * 


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
d = r2C - r1A
fs = pp_s(d, k, L0)
fsa, fsc = fs*uy, -fs*uy

ecchime = 1 