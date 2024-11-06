# inUvectors
"""
Created on Thu Sep 28 15:18:09 2023

@author: simone.lucertini
"""

import Unit_struct as Us


U1 = dict(Us.Unit)
U1["Bindex"] = 0
U1["ulocal"]  = [ 1.0, 0]

U2 = dict(Us.Unit)
U2["Bindex"]= 1
U2["ulocal"]= [ 1.0, 0]

Uvectors = [U1, U2]
del (U1, U2)
