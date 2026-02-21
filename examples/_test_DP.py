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
double_pendulum = PlanarMultibodyModel()
T, uT = double_pendulum.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

# #%% Save results
# import os
# output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', '_test_DP.txt')
# nB = double_pendulum.nB
# nC = double_pendulum.nC
# nB3 = nB * 3
# header = '\t'.join(['t'] + [f'B{i+1}_{c}' for i in range(nB) for c in ['x','y','p']])
# np.savetxt(output_file, np.column_stack([T, uT[:, :nB3]]),
#            delimiter='\t', header=header, comments='', fmt='%.8f')
# print(f"[_test_DP] Done. nB={nB}, nC={nC}, DOF={nB*3-nC}, points={len(T)}")
# print(f"  Results: {output_file}")

plot_comparison(T, uT, matlab_filename='DP.txt', model_title='DP')
