# body properties definition
import numpy as np

#def inBodies():
import Body_struct as Bs


# build Body 1 structure
B1 = dict(Bs.Body) #Copy existing (empty" DICT)
B1["m"] = 5.0
B1["J"] = 4.0
B1["r"] = [1.0, 0.2]

# build Body 1 structure
B2 = dict(Bs.Body) #Copy existing (empty" DICT)
B2["m"] = 2.0
B2["J"] = 0.2
B2["r"] = [1.25, -0.233]
B2["p"] = np.pi/6

# build full body structure
Bodies = [B1, B2] # LIST of DICTS
del (B1, B2)

with open("inAnimate.py") as f: # !!! originariament era in inData
    exec(f.read())   

    