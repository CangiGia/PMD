import numpy as np
import scipy as sc
import matplotlib as mpl
from PMD.src.builder import *
from PMD.src.functions import *
from PMD.src.solver import *
import matplotlib.pyplot as plt


#%% bodies
b1 = Body(m=1, J=0.0833, r=[0.0,0.5], p=0.0)
b2 = Body(m=2, J=0.375, r=[0.75,1.0], p=0.0)

#%% points
p0 = Point(Bindex=0, sPlocal=np.array([0.0, 0.0]))     #// revolute joint between b1 and b0 - ground side (b0 is always the ground)
p1 = Point(Bindex=1, sPlocal=np.array([0.0, -0.5]))    #// revolute joint between b1 and b0 - body 1 side (b0 is always the ground)
p2 = Point(Bindex=1, sPlocal=np.array([0.0, 0.5]))     #// revolute joint between b1 and b2 - body 1 side
p3 = Point(Bindex=2, sPlocal=np.array([-0.75, 0]))     #// revolute joint between b1 and b2 - body 2 side

#%% joints
j1 = Joint(type="rev", iPindex=0, jPindex=1) #// Revolute joint between b1 and b0
j2 = Joint(type="rev", iPindex=2, jPindex=3) #// Revolute joint between b1 and b2

#%% forces
s3 = Force()
s3 = Force(type="weight") #// only weight force, acting along -y axis

#%% double pendulum model creation
double_pendulum = PlanarDynamicModel()
time, solution = double_pendulum.solve(method='Radau')

fig, ax = plt.subplots(2, 1, figsize=(12, 12), sharex=True)

ax[0].plot(time, solution[:, 0], label=r'$x_1$', color='b', linestyle='-', linewidth=2)
ax[0].plot(time, solution[:, 1], label=r'$y_1$', color='r', linestyle='-', linewidth=2)
ax[0].set_ylabel('Displacement (m)', fontsize=12)
ax[0].legend()

ax[1].plot(time, solution[:, 3], label=r'$x_2$', color='b', linestyle='-', linewidth=2)
ax[1].plot(time, solution[:, 4], label=r'$y_2$', color='r', linestyle='-', linewidth=2)
ax[1].set_ylabel('Displacement (m)', fontsize=12)
ax[1].legend()

# ax[2].plot(time, solution[:, 2], label=r'$\psi_1$', color='b', linestyle='-', linewidth=2)
# ax[2].plot(time, solution[:, 5], label=r'$\psi_2$', color='r', linestyle='-', linewidth=2)
# ax[2].set_xlabel('Time (s)', fontsize=12)
# ax[2].set_ylabel('Displacement (rad)', fontsize=12)
# ax[2].legend()

fig.suptitle('System responses', fontsize=14)
plt.show()