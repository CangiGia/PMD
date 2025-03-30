import numpy as np
import scipy as sc
from PMD.src.builder import *
from PMD.src.functions import *
from PMD.src.solver import *
import matplotlib.pyplot as plt


#%% bodies
b1 = Body(m=8, J=2e-03, r=[0.25,-0.25], p=(np.pi/4))
b2 = Body(m=2, J=5e-04, r=[0.525,-0.425], p=(-(70.404)*np.pi/180))

#%% points
p0 = Point(Bindex=0, sPlocal=np.array([0, 0]))         # revolute joint between b1 and b0 - ground side
p1 = Point(Bindex=1, sPlocal=np.array([-0.25, 0]))     # revolute joint between b1 and b0 - body 1 side
p2 = Point(Bindex=1, sPlocal=np.array([0.25, 0]))      # revolute joint between b1 and b2 - body 1 side
p3 = Point(Bindex=2, sPlocal=np.array([-0.079, 0]))    # revolute joint between b1 and b2 - body 2 side

#%% joints
j1 = Joint() #// Revolute joint between b1 and b0
j1.type = 'rev'
j1.iPindex = 0
j1.jPindex = 1

j2 = Joint() #// Revolute joint between b1 and b0
j2.type = 'rev'
j2.iPindex = 2
j2.jPindex = 3

#%% forces
s3 = Force()
s3.type = 'weight'
s3.wgt = np.array([0, -1])

#%% double pendulum model creation
my_dynamic_model = PlanarDynamicModel(verbose=False)
time, solution = my_dynamic_model.solve(method='RK45')

plt.figure()
plt.plot(time, solution[:,3])
plt.show()
