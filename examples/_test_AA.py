import numpy as np
import scipy as sc
from PMD.src.builder import *
from PMD.src.solver import *
from PMD.src.mechanics import *
import matplotlib.pyplot as plt

#%% bodies
B1 = Body(m=2, J=0.5, r=[0.4398, 0.2512], p=6.2465)   # lower suspension arm
B2 = Body(m=30, J=2.5, r=[0.6817, 0.3498], p=0.0783)  # wheel assembly
B3 = Body(m=1, J=0.5, r=[0.4463, 0.4308], p=6.5222)   # upper suspension arm

#%% points
q1 = Point(Bindex=1, sPlocal=np.array([-0.24, 0.0]))   # Q point - lower suspension arm side
a1 = Point(Bindex=1, sPlocal=np.array([0.18, 0.0]))    # A point - lower suspension arm side
a2 = Point(Bindex=2, sPlocal=np.array([-0.07, -0.10])) # A point - wheel assembly side
b2 = Point(Bindex=2, sPlocal=np.array([-0.10, 0.12]))  # B point - wheel assembly side
b3 = Point(Bindex=3, sPlocal=np.array([0.13, 0.0]))    # B point - upper suspension arm side
o3 = Point(Bindex=3, sPlocal=np.array([-0.13, 0.0]))   # O point - upper suspension arm side
o0 = Point(Bindex=0, sPlocal=np.array([0.32, 0.40]))   # O point - ground
q0 = Point(Bindex=0, sPlocal=np.array([0.20, 0.26]))   # Q point - ground
e1 = Point(Bindex=1, sPlocal=np.array([0.0, 0.0]))     # E point - Lower suspension arm side
f0 = Point(Bindex=0, sPlocal=np.array([0.38, 0.43]))   # F point - Ground

#%% joints
j1 = Joint(type='rev', iPindex=0, jPindex=7)  # Revolute joint in Q
j2 = Joint(type='rev', iPindex=1, jPindex=2)  # Revolute joint in A
j3 = Joint(type='rev', iPindex=3, jPindex=4)  # Revolute joint in B
j4 = Joint(type='rev', iPindex=5, jPindex=6)  # Revolute joint in O

#%% forces
s1 = Force(type='ptp', iPindex=8, jPindex=9, k=90000, L0=0.23, dc=1100)  # spring-damper force

# custom wheel contact force
def my_force(B2, s2):
    """Custom force used to define the wheel contact condition."""
    dely = B2.r[1] - s2.L0
    if dely < 0:
        fy = (s2.k * dely + s2.dc * B2.dr[1]).item()
        fsd = np.array([0, -fy])
        B2._f = B2._f + fsd.reshape(2,1)

s2 = Force(type='user', k=50000, L0=0.35, dc=1000)  # wheel contact force
s2.callback = my_force

s3 = Force(type='weight')  # gravity force

#%% solution
quarter_car = PlanarMultibodyModel()
T, uT = quarter_car.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001),
                          ic_correct=True)

#%% Save results
import os
output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', '_test_AA.txt')
nB = quarter_car.nB
nC = quarter_car.nC
nB3 = nB * 3
header = '\t'.join(['t'] + [f'B{i+1}_{c}' for i in range(nB) for c in ['x','y','p']])
np.savetxt(output_file, np.column_stack([T, uT[:, :nB3]]),
           delimiter='\t', header=header, comments='', fmt='%.8f')
print(f"[_test_AA] Done. nB={nB}, nC={nC}, DOF={nB*3-nC}, points={len(T)}")
print(f"  Results: {output_file}")
