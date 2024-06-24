import sys
import os
import numpy as np

sys.path.append(os.path.abspath("."))
from useful_code.numerov_method import *

def f(E):
    return lambda r: -2*(E + 1/r)

def uAsym(r):
    return r*np.exp(-r) 

def getRMin(E):
    rs,us = runNumerov(rStart,h,uAsym(rStart),uAsym(rStart+h),f(E),N)
    return us[-1]

rStart = 100
h = -0.1
N = 999

ERange = -np.arange(0.01,2,0.01)
print(min(ERange,key=getRMin))
