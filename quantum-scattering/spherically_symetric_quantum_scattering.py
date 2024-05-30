import sys
import os
import numpy as np
from scipy import constants

sys.path.append(os.path.abspath("."))
from useful_code.numerov_method import *

class ScatteringSystem:
    def __init__(self,E,l,m,V,r_0,h,u_0,u_1,r_end,simpleUnits=False):
        #system parameters
        self.E = E
        self.l = l
        self.m = m
        self.V = V
        self.simpleUnits = simpleUnits

        #simulation parameters
        self.r_0 = r_0
        self.h = h
        self.N = int(np.ceil((r_end-r_0)/h))

        self.rList,self.uList = runNumerov(r_0,h,u_0,u_1,self.F,self.N)

    def F(self,r):
        #removes certain divergences that don't actually occur
        if r == 0 and self.r_0 == 0:
            return 0
        
        if self.simpleUnits:
            return self.V(r) + self.l*(self.l+1)/(r**2) - self.E
        
        return 2*self.m*self.V(r)/(constants.hbar**2) + self.l*(self.l+1)/(r**2) - 2*self.m*self.E/(constants.hbar**2)


    def findNearestu(self,r):
        index = np.floor(r/self.h)
        return self.uList[index]

        
        


