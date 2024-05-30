import sys
import os
from scipy import constants
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
r_0 = 0.5*rho
r_1 = 0.8*rho
h = r_1 - r_0
u_0 = uSmallr(r_0)
u_1 = uSmallr(r_1)

#system parameters
E = 1
r_end = 10000

#run simulation
scatteringSys = ScatteringSystem(E,V_LJ,r_0,h,u_0,u_1,r_end,6,simpleUnits=False)

print(scatteringSys.totalCrossSection)