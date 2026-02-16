import numpy as np
import matplotlib.pyplot as plt
from utils import *

## Indata (Geometry, meterial)
E = 210e9
L = 2
A0 = 78.5e-4 #cm3
A_spec = A0*2
P = 150e3
sigmas = 230e6

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


Ex, Ey = coordxtr(Edof, Coord, Dofs)    # x och y-koordinater för varje element, alltså start och slutpunkt för aktuell stång. 

#Hjälpvariabler:
nel = Edof.shape[0]  # Antal element
ndofs = 14 #Totalt antal frihetsgrader

# Fördefinera styvhetsmatrisen och kraftvektorn
K = np.zeros((ndofs, ndofs))  #alltså är K styvhetsmatrisen som är 14*14 och f är kraftvektorn som är 14 lång. 
f = np.zeros(ndofs)

# Assemblera elemented
for el in range(nel):
    x1, x2 = Ex[el]
    y1, y2 = Ey[el]

    L_el = np.sqrt((x2-x1)**2 + (y2-y1)**2) #här beräknar vi hur långt elementet (stången) är mellan både y och x coordinaterna för noderna. 
    c = (x2-x1) / L_el #cossinus, det bärknas genom att ställa upp rätvinkliga triangeln och ta formeln för cos och sinus. Det är alltså därifrån denna definering kommer ifrån. 
    s = (y2-y1) / L_el  #tex sin blir skillnaden i y koordinaterna (motstående) / längden av det aktuella stången (hypotenusa. )
    
    A = A_spec if el in [3, 5, 8] else A0 #alltså tvärsnittsarean blir a_spec om det är en av de stänger som ska ha en annan tvärsnittsarea. 

    k_faktor = (E*A/L_el)  # detta är som formeln vi gick igenom på föreläsningen, alltså aktuell tväsnittarea och aktuellt längd på elementet (stålpen) vi just nu kollar på

    Ke = k_faktor*np.array([  #detta är Ekvation 11.18 i kursboken och det som tagitsfram här är produkten med k faktorn och stychetsmatrisen
        [c**2, c*s, -c**2, -c*s],
        [c*s, s**2, -c*s, -s**2],
        [-c**2, -c*s, c**2, c*s],
        [-c*s, -s**2, c*s, s**2]
    ])
    
    #Assemblera in element styvhetsmatrisen och globala matrisen
    K = assem(Edof[el, :], K, Ke)

# Lägg till kraften P i lastvektorn:
f[11] = -P

bcdofs = np.array([1, 2, 3, 4]) #detta är frihetsgraderna för nod 1 och 2, vi vill göra detta för att det berättar om det finns några frihetsgrader som är låsta
bcvals = np.array([0, 0, 0, 0]) #denna delen talar om hur mycket de låsta frihetsgraderna får röra på sig, och eftersom dessa är fasta stöd så rör de sig inte alls, och vi sätter värdena till 0

# Lös ekvations systemet (: använd solveq i utils.py)
a, Q = solveq(K, f, bcdofs, bcvals) 


# Plotta deformerad mesh (: använd eldisp2 i utils.py)
print("Knutförskjutningar (a):")
print(a) #printar förkjutningsvektorn. 

plt.figure()
eldraw2(Ex, Ey, width=1, color="black") # Ursprunglig form
Ed = extract_eldisp(Edof, a)
eldisp2(Ex, Ey, Ed, sfac=50, width=2, color="blue") # Deformerad form

plt.title("Fackverk: Odefomerad (svart) och Deformerad (blå, 50x skala)")
plt.axis("equal")
plt.show()

# Räkna ut krafter och spänningar i varje element
Stångkrafter = []
Spänningar = []

for el in range(nel):
    #hämta ut det aktuella elementets (stångens) koordinater, 
    ex = Ex[el]
    ey = Ey[el]
    #formel för spänningar är sigma = N/A, där N är stångkraften och A är tvärsnittsarean. 
    
    A = A_spec if el in [3, 5, 8] else A0 # ser till så att vi får rätt are för stången vi är på nu. 
    ep = [E, A] #detta står i bar2s funktionen i utils. 
     #vi hämtar förskutningarna för varje element från den globala vektorn a genom att använda extract_eldisp ifrån utils.
     #edof[el,:] bestämmet vilka frihetsgrader denna stången har. 
    ed = extract_eldisp(Edof[el:el+1, :], a)

    #anbvänd bars2s funktion för att beräkna spänningen i aktuellt element. 
    N_el = bar2s(ex, ey, ep, ed) #denna returnerar en array som innehåller stångkrafterna. 

    #stångkraften beräknas sedan enligt formeln som jag skrev tidigare i dokumentet, 
    sigma_el =  N_el/ A 
    Stångkrafter.append(N_el[0])  #lägga till stångkrafter och spänningar i tidigare listor
    Spänningar.append(sigma_el[0])

    #print(Stångkrafter, Spänningar)
    print(f"{el+1} {N_el[0]/1000:10.2f}   {sigma_el[0]/1e6:10.2f}") #printa ut alla sigam och stångkrafter samt avrunda till finare tal. 


Stångkrafter = np.array(Stångkrafter)
Spänningar = np.array(Spänningar)

imax = np.argmax(Spänningar)   # störst drag
imin = np.argmin(Spänningar)   # störst tryck

print("\nResultat:")
print(f"Störst dragspänning:  Stång {imax+1}  sigma = {Spänningar[imax]/1e6:.2f} MPa")
print(f"Störst tryckspänning: Stång {imin+1}  sigma = {Spänningar[imin]/1e6:.2f} MPa")


    #... tips: använd bar2s i utils.py
