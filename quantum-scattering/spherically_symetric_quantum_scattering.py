import sys
import os
import numpy as np
from scipy import constants
from scipy import special

sys.path.append(os.path.abspath("."))
from useful_code.numerov_method import *

def findPhaseShift(r_max,u_max,r_end,u_end,l,k):
    K = (r_max*u_end)/(r_end*u_max)
    numerator = K*special.spherical_jn(l,k*r_max) - special.spherical_jn(l,k*r_end)
    denominator = K*special.spherical_yn(l,k*r_max) - special.spherical_yn(l,k*r_end)
    return np.arctan(numerator/denominator)


class ScatteringSystem:
    def __init__(self,E,V,r_0,h,u_0,u_1,r_max,lCutoff,simpleUnits=False):
        #system parameters
        self.E = E
        self.k = np.sqrt(6.12*E)
        self.lambbda = 2*np.pi/self.k
        self.V = V
        self.simpleUnits = simpleUnits

        #cutoff for sums in l
        self.lCutOff = lCutoff

        #simulation parameters
        self.h = h
        self.N = int(np.ceil((r_max + self.lambbda - r_0)/h))

        #intial and final r values
        self.r_0 = r_0
        self.r_maxIndex = int(np.ceil((r_max-r_0)/h))
        self.r_max = h*self.r_maxIndex + r_0
        self.r_end = r_0 + self.N*h

        self.uLists = []
        for l in range(self.lCutOff):
            F = lambda r: self.FL(r,l)
            self.rList, ulList = runNumerov(r_0,h,u_0,u_1,F,self.N)
            self.uLists.append(ulList)

        self.findPhaseShifts()
        self.findTotalCrossSection()

    def FL(self,r,l):
        #removes certain divergences that don't actually occur
        if r == 0 and self.r_0 == 0:
            return 0
        
        if self.simpleUnits:
            return self.V(r) + l*(l+1)/(r**2) - self.E
        
        #choose units of p^2 for E and V is 
        return 6.12*self.V(r) + l*(l+1)/(r**2) - 6.12*self.E

    def findPhaseShifts(self):
        self.phaseShifts = []
        for l in range(self.lCutOff):
            u_max = self.uLists[l][self.r_maxIndex]
            u_end = self.uLists[l][-1]
            self.phaseShifts.append(findPhaseShift(self.r_max,u_max,
                                           self.r_end,u_end,l,self.k))
            
    def findTotalCrossSection(self):
        elementsOfSum = [(2*l + 1)*(np.sin(deltaL)**2) for l,deltaL in enumerate(self.phaseShifts)]
        
        self.totalCrossSection = 4*np.pi*sum(elementsOfSum)/(self.k**2)











        

        
        


