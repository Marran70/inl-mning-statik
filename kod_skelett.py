import numpy as np
import matplotlib.pyplot as plt
from utils import *

## Indata (Geometry, meterial)
E = ...
A0 = ...
P = ...
...

## Topology
Edof = np.array([
    [dofx1, dofy1, dofx2, dofy2],
    ...
], dtype=int )

# Koordinater för varje nod:
Coord = np.array([
    [0.0, 0.0],
    ...
])

# x-koordinater för varje element
Ex = np.array([
    [0.0, 0.0],
    ...
])

# y-koordinater för varje element
Ey = np.array([
    [0.0, 0.0],
    ...
])

#Hjälpvariabler:
nel = ...  # Antal element
ndofs = ... #Totalt antal frihetsgrader

# Plot mesh (tips: använd eldraw2 i utils.py)
...

# Fördefinera styvhetsmatrisen och kraftvektorn
K = ...
f = ...

# Assemblera elemented
for el in range(nel):
    
    #Räkna ut styvhetsmatrisen (Se föreläsningar eller Ekvation 11.18 i kursboken [Hållfasthetslära, Allmänna tillstånd])
    Ke = ...
    
    #Assemblera in element styvhetsmatrisen och globala matrisen
    K = assem(Edof[el, :], K, Ke)

# Lägg till kraften P i lastvektorn:
...


# Lös ekvations systemet (: använd solveq i utils.py)
...


# Plotta deformerad mesh (: använd eldisp2 i utils.py)
...

# Räkna ut krafter och spänningar i varje element
for el in range(nel):
    #... tips: använd bar2s i utils.py

#osv...
