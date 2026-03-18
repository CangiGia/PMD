import numpy as np
import scipy as sc
from PMD.src import *
import matplotlib.pyplot as plt
from PMD.examples._plot_utils import plot_comparison

#%% bodies
B1 = Body(m=2, J=0.5, r=[0.4398, 0.2512], p=-0.0367)   # lower suspension arm
B2 = Body(m=30, J=2.5, r=[0.6817, 0.3498], p=0.0783)  # wheel assembly
B3 = Body(m=1, J=0.5, r=[0.4463, 0.4308], p=6.5222)   # upper suspension arm

#%% markers
q1 = B1.add_marker([-0.24, 0.0])       # Q point - lower suspension arm side
a1 = B1.add_marker([0.18, 0.0])        # A point - lower suspension arm side
a2 = B2.add_marker([-0.07, -0.10])     # A point - wheel assembly side
b2 = B2.add_marker([-0.10, 0.12])      # B point - wheel assembly side
b3 = B3.add_marker([0.13, 0.0])        # B point - upper suspension arm side
o3 = B3.add_marker([-0.13, 0.0])       # O point - upper suspension arm side
o0 = Ground.add_marker([0.32, 0.40])   # O point - ground
q0 = Ground.add_marker([0.20, 0.26])   # Q point - ground
e1 = B1.add_marker([0.0, 0.0])         # E point - Lower suspension arm side
f0 = Ground.add_marker([0.38, 0.43])   # F point - Ground

#%% joints
j1 = Joint(type='rev', iMarker=q1, jMarker=q0)  # Revolute joint in Q
j2 = Joint(type='rev', iMarker=a1, jMarker=a2)  # Revolute joint in A
j3 = Joint(type='rev', iMarker=b2, jMarker=b3)  # Revolute joint in B
j4 = Joint(type='rev', iMarker=o3, jMarker=o0)  # Revolute joint in O

#%% forces
s1 = Force(type='ptp', iMarker=e1, jMarker=f0, k=90000, L0=0.23, dc=1100)  # spring-damper force

# custom wheel contact force
s2 = Force(type='user', k=50000, L0=0.35, dc=1000)  # wheel contact force

def my_force():
    """Custom force used to define the wheel contact condition."""
    dely = B2.r[1] - s2.L0
    if dely < 0:
        fy = (s2.k * dely + s2.dc * B2.dr[1]).item()
        fsd = np.array([0, -fy])
        B2._f = B2._f + fsd.reshape(2,1)

s2.callback = my_force

s3 = Force(type='weight')  # gravity force

#%% solution
quarter_car = PlanarMultibodyModel(
    bodies=[B1, B2, B3],
    joints=[j1, j2, j3, j4],
    forces=[s1, s2, s3])
T, uT = quarter_car.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001),
                          ic_correct=True)

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='AA.txt', model_title='AA')
