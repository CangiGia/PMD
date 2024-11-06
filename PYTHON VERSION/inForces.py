# inForces
"""
Created on Thu Sep 28 15:53:51 2023

@author: simone.lucertini
"""

import Force_struct as Fs

F1 = dict(Fs.Force)
F1["type"] = "ptp"; # default
F1["iPindex"] = 2
F1["jPindex"] = 1
F1["k"] = 20
F1["L0"] = 0.6

F2 = dict(Fs.Force)
F2["type"] = "weight"  # include the weight
F2["gravity"]= 9.81   #default
F2["wgt"] = [0, -1] # default

Forces = [F1, F2]
del (F1, F2)


 