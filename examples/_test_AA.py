import numpy as np
import scipy as sc
from PMD.src.builder import *
from PMD.src.solver import *
from PMD.src.functions import *
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

#%% solution and plotting
quarter_car = PlanarMultibodyModel()
time, solution = quarter_car.solve(method='Radau')
# reactions = quarter_car.get_reactions()

plt.figure(figsize=(15, 5))
plt.plot(time, solution[:,4], label='Wheel vertical position', color='y', linestyle='-', linewidth=2)
plt.ylabel('Displacement (m)', fontsize=12)
plt.legend()
plt.title('Wheel Assembly Response', fontsize=14)
plt.show()

# plt.figure(figsize=(15, 5))
# plt.plot(time, reactions[:, 0], label='Reaction force at Q', color='b', linestyle='-', linewidth=2)
# plt.ylabel('Force (N)', fontsize=12)
# plt.legend()
# plt.title('Reaction force at Q', fontsize=14)
# plt.show()