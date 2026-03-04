import numpy as np
import scipy as sc
import matplotlib as mpl
from PMD.src.builder import *
from PMD.src.mechanics import *
from PMD.src.solver import *
import matplotlib.pyplot as plt
plt.style.use('dark_background')
from PMD.examples._plot_utils import plot_comparison


#%% bodies
b1 = Body(m=1, J=0.0833, r=[0.0,0.5], p=0.0)
b2 = Body(m=2, J=0.375, r=[0.75,1.0], p=0.0)

#%% points
p0 = Point(body=Ground, sPlocal=np.array([0.0, 0.0]))     #// revolute joint between b1 and b0
p1 = Point(body=b1, sPlocal=np.array([0.0, -0.5]))    #// revolute joint between b1 and b0
p2 = Point(body=b1, sPlocal=np.array([0.0, 0.5]))     #// revolute joint between b1 and b2
p3 = Point(body=b2, sPlocal=np.array([-0.75, 0]))     #// revolute joint between b1 and b2

#%% joints
j1 = Joint(type="rev", iPoint=p0, jPoint=p1) #// Revolute joint between b1 and b0
j2 = Joint(type="rev", iPoint=p2, jPoint=p3) #// Revolute joint between b1 and b2

#%% forces
s3 = Force()
s3 = Force(type="weight") #// only weight force, acting along -y axis

#%% double pendulum model creation
double_pendulum = PlanarMultibodyModel(
    bodies=[b1, b2],
    joints=[j1, j2],
    forces=[s3],
    points=[p0, p1, p2, p3])
T, uT = double_pendulum.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='DP.txt', model_title='DP')
