# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 13:50:59 2023

@author: simone.lucertini
"""

import Point_struct as Ps

P1 = dict(Ps.Point)
P1["Bindex"] = 0
P1["sPlocal"] = [ 0, 0.2]

P2 = dict(Ps.Point)
P2["Bindex"] = 1
P2["sPlocal"] = [ 0, 0]

P3 = dict(Ps.Point)
P3["Bindex"] = 2
P3["sPlocal"] = [ 0, 0.5]

Points = [P1, P2, P3]
del (P1, P2, P3)