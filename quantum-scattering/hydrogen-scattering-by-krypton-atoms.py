import sys
import os
import numpy
from scipy import constants
# eg. constants.hbar

sys.path.append(os.path.abspath("."))

from useful_code.numerov_method import *

def V(r):
    return 0

def F(l,r,E):
    #if hbar^2/2m != 1 then units need to change here
    return V(r) + (constants.hbar**2)*l*(l+1)
