import numpy as np
import scipy as sc
from PMD.src.builder import *
from PMD.src.functions import *
from PMD.src.solver import *
import matplotlib.pyplot as plt


#%% bodies
b1 = Body(m=4.2018301966, J=9.9897988328E-04, r=[0.25,0.0], p=0.0)
b2 = Body(m=1.0814301966, J=2.4488321698E-04, r=[0.525,0.0], p=0.0)

#%% points
p0 = Point(Bindex=0, sPlocal=np.array([0, 0]))         # revolute joint between b1 and b0 - ground side
p1 = Point(Bindex=1, sPlocal=np.array([-0.25, 0]))     # revolute joint between b1 and b0 - body 1 side
p2 = Point(Bindex=1, sPlocal=np.array([0.25, 0]))      # revolute joint between b1 and b2 - body 1 side
p3 = Point(Bindex=2, sPlocal=np.array([-0.05, 0]))     # revolute joint between b1 and b2 - body 2 side

#%% joints
j1 = Joint() #// Revolute joint between b1 and b0
j1.type = 'rev'
j1.iPindex = 0
j1.jPindex = 1

j2 = Joint() #// Revolute joint between b1 and b2
j2.type = 'rev'
j2.iPindex = 2
j2.jPindex = 3

#%% forces
s3 = Force()
s3.type = 'weight'
s3.wgt = np.array([0, -1])

#%% double pendulum model creation
my_dynamic_model = PlanarDynamicModel(verbose=False)
time, solution = my_dynamic_model.solve(method='BDF')

plt.figure()
plt.plot(time, solution[:,0])
plt.show()

np.savetxt('uT_python.txt', solution)
np.savetxt('T_python.txt', time)
