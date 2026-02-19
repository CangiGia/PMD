import numpy as np
import scipy as sc
from PMD.src.builder import *
from PMD.src.mechanics import *
from PMD.src.solver import *
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

#%% model simulation
my_dynamic_model = PlanarMultibodyModel()
T, uT = my_dynamic_model.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

#%% Save results
import os
output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', '_test_SP.txt')
nB = my_dynamic_model.nB
nC = my_dynamic_model.nC
nB3 = nB * 3
header = '\t'.join(['t'] + [f'B{i+1}_{c}' for i in range(nB) for c in ['x','y','p']])
np.savetxt(output_file, np.column_stack([T, uT[:, :nB3]]),
           delimiter='\t', header=header, comments='', fmt='%.8f')
print(f"[_test_SP] Done. nB={nB}, nC={nC}, DOF={nB*3-nC}, points={len(T)}")
print(f"  Results: {output_file}")
