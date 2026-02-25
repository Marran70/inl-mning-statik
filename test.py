import numpy as np
import matplotlib.pyplot as plt
from utils import *

## Indata (Geometry, meterial)
E = 210e9
L = 2
A0 = 78.5e-4 
A_spec = A0*2
P = 150e3
sigma_s = 230e6

## Topology
Edof = np.array([  #det är denna som bestämmer hur figuren ser ut. den fungerar så att det är frihetsgraderna för startnoden samt slutnoden.
    [3, 4, 7, 8],     # Stång 1 
    [3, 4, 5, 6],     # Stång 2 
    [1, 2, 7, 8],     # Stång 3 
    [1, 2, 5, 6],     # Stång 4 (Area: 2A0) 
    [5, 6, 7, 8],     # Stång 5 
    [5, 6, 9, 10],    # Stång 6 (Area: 2A0) 
    [7, 8, 9, 10],    # Stång 7 
    [7, 8, 11, 12],   # Stång 8 
    [9, 10, 11, 12],  # Stång 9 (Area: 2A0) 
    [9, 10, 13, 14],  # Stång 10 
    [11, 12, 13, 14]  # stång 11
], dtype=int )


Coord = np.array([    # Koordinater för varje nod:
    [0.0, 0.0],           # Nod 1 (Origo) 
    [L, 0.0],             # Nod 2 
    [0.0, 2*L],           # Nod 3 
    [L, 2*L],             # Nod 4 
    [L, 3*L],             # Nod 5 
    [3*L, 2*L],           # Nod 6 
    [3*L, 3*L]            # Nod 7 
])

Dofs = np.array([ #frihetsgrader för varje nod. 
    [1, 2], #nod 1
    [3, 4],  #nod 2 osv...
    [5, 6], 
    [7, 8], 
    [9, 10], 
    [11, 12], 
    [13, 14]
])


Ex, Ey = coordxtr(Edof, Coord, Dofs) 

print(Ex)
print(Ey)