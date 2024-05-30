import sys
import os
from scipy import constants
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import numpy as np
from spherically_symetric_quantum_scattering import ScatteringSystem

# some constants: energy is measure in meV and distance in units of rho
# means 2m/hbar^2 = 6.12 meV^-1 rho^-2
epsilon = 5.9
rho = 1
alpha = 6.12
C = np.sqrt(epsilon*alpha/25)

def V_LJ(r):
    return epsilon*((rho/r)**12 - 2*((rho/r)**6))

def uSmallr(r):
    return np.exp(-C/(r**5))


#find intial conditions and h
r_0 = 0.7*rho
r_1 = 0.701*rho
h = r_1 - r_0
u_0 = uSmallr(r_0)
u_1 = uSmallr(r_1)

#system parameters
EArray = np.arange(.1,3.5,.01)
r_end = 5

#run simulation
lst = []
for E in EArray:
    lst.append(ScatteringSystem(E,V_LJ,r_0,h,u_0,u_1,r_end,8,simpleUnits=False).totalCrossSection)
    print(E)

#filterlst = savgol_filter(lst, 3, 3) 

plt.plot(EArray,lst)
#plt.plot(EArray,filterlst)

plt.show()

#peak at .5 1.5 and 2.5