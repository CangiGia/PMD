# inAnimate
"""
Created on Fri Sep 29 09:57:19 2023

@author: simone.lucertini
"""
# C2 = Point_struct;
# C2.Bindex = 2;
# C2.sPlocal = [0;-0.5];

# A1 = Point_struct;
# A1.Bindex = 1;
# A1.sPlocal = [0.2;0];
# 
# B1 = Point_struct;
# B1.Bindex = 1;
# B1.sPlocal = [-0.2;0];

from inBodies import Bodies

Points_anim = []

Bodies[0]["shape"] = "rect"
Bodies[0]["W"] = 0.4
Bodies[0]["H"] = 0.2

Bodies[1]["shape"]= "line"
Bodies[1]["H"] = 1.0


# Parameters for defining animation/plot axes 
xmin = -0.5
xmax =  2.0
ymin = -1.0
ymax =  0.5


