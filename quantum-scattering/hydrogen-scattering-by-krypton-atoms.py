import sys
import os
from scipy import constants
from spherically_symetric_quantum_scattering import ScatteringSystem

# some constants
epsilon = 5.9*(10**-3)*constants.elementary_charge
rho = 3.67*constants.angstrom

def V_LJ(r):
    return epsilon*((rho/r)**12 - 2*((rho/r)**6))

