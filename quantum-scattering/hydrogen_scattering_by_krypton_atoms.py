import sys
import os
import numpy
from scipy import constants
# eg. constants.hbar

sys.path.append(os.path.abspath("."))

from useful_code.numerov_method import *

# some constants
epsilon = 5.9*(10**-3)*constants.elementary_charge
rho = 3.67*constants.angstrom

def V_LJ(r):
    return epsilon*((rho/r)**12 - 2*((rho/r)**6))


def F(l,r,E,m,V,units=False):
    if r == 0:
        #allowed to do this as it only occurs at u_0 = 0, so it has no effect (make this clearer)
        return 0
    if units:
        return V(r) + l*(l+1)/(r**2) - E
    return 2*m*V(r)/(constants.hbar**2) + l*(l+1)/(r**2) - 2*m*E/(constants.hbar**2)





