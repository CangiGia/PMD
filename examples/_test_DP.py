import numpy as np
import scipy as sc
import matplotlib as mpl
from PMD.src import *
import matplotlib.pyplot as plt
plt.style.use('dark_background')
from PMD.examples._plot_utils import plot_comparison


#%% bodies
b1 = Body(m=1, J=0.0833, r=[0.0,0.5], p=0.0)
b2 = Body(m=2, J=0.375, r=[0.75,1.0], p=0.0)

#%% markers
p0 = Ground.add_marker([0.0, 0.0])    #// revolute joint between b1 and b0
p1 = b1.add_marker([0.0, -0.5])       #// revolute joint between b1 and b0
p2 = b1.add_marker([0.0, 0.5])        #// revolute joint between b1 and b2
p3 = b2.add_marker([-0.75, 0])        #// revolute joint between b1 and b2

#%% joints
j1 = Joint(type="rev", iMarker=p0, jMarker=p1) #// Revolute joint between b1 and b0
j2 = Joint(type="rev", iMarker=p2, jMarker=p3) #// Revolute joint between b1 and b2

#%% forces
s3 = Force(type="weight") #// only weight force, acting along -y axis

#%% double pendulum model creation
double_pendulum = PlanarMultibodyModel(
    bodies=[b1, b2],
    joints=[j1, j2],
    forces=[s3])
T, uT = double_pendulum.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='DP.txt', model_title='DP')
