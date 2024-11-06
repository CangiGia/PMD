# inJoints
"""
Created on Thu Sep 28 16:11:30 2023

@author: simone.lucertini
"""


import Joint_struct as Js

J1 = dict(Js.Joint)
J1["type"] = "tran"
J1["iPindex"]= 2
J1["jPindex"]= 1
J1["iUindex"] = 2
J1["jUindex"]= 1

J2 = dict(Js.Joint)
J2["type"] = "rev"
J2["iPindex"] = 2
J2["jPindex"] = 3

Joints = [J1, J2]
del (J1, J2) #delete unnecessary vars.
