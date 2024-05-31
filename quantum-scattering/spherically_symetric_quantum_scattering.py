import sys
import os
import numpy as np
from scipy import special

sys.path.append(os.path.abspath("."))
from useful_code.numerov_method import *

#units, defines 2m/hbar^2 units
units = 6.12

class ScatteringSystem:
    def __init__(self,E,V,r_0,h,u_0,u_1,r_max,lCutoff):
        #system parameters
        self.E = E
        self.k = np.sqrt(units*E)
        self.lambbda = 2*np.pi/self.k
        self.V = V

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

        #runs simulations
        self.uLists = []
        for l in range(self.lCutOff):
            F = lambda r: self.FL(r,l)
            self.rList, ulList = runNumerov(r_0,h,u_0,u_1,F,self.N)
            self.uLists.append(ulList)

        #finds phase shifts and cross sections
        self.findPhaseShifts()
        self.findTotalCrossSection()

    def FL(self,r,l):
        #removes certain divergences that don't actually occur
        if r == 0 and self.r_0 == 0:
            return 0
        return units*self.V(r) + l*(l+1)/(r**2) - units*self.E
    
    def findlPhaseShift(self,l,u_max,u_end):
        K = (self.r_max*u_end)/(self.r_end*u_max)
        numerator = K*special.spherical_jn(l,self.k*self.r_max) - special.spherical_jn(l,self.k*self.r_end)
        denominator = K*special.spherical_yn(l,self.k*self.r_max) - special.spherical_yn(l,self.k*self.r_end)
        return np.arctan(numerator/denominator)


    def findPhaseShifts(self):
        self.phaseShifts = []
        for l in range(self.lCutOff):
            u_max = self.uLists[l][self.r_maxIndex]
            u_end = self.uLists[l][-1]
            self.phaseShifts.append(self.findlPhaseShift(l,u_max,u_end))
            
    def findTotalCrossSection(self):
        elementsOfSum = [(2*l + 1)*(np.sin(deltaL)**2) for l,deltaL in enumerate(self.phaseShifts)]
        self.totalCrossSection = 4*np.pi*sum(elementsOfSum)/(self.k**2)

    def differentialCrossSection(self,theta):
        termOfSum = lambda deltaL,l:  np.abs((2*l+1)*np.exp(complex(0,1)*deltaL)*np.sin(deltaL)*
                                special.eval_legendre(l,np.cos(theta)))**2
        elementsOfSum = [termOfSum(deltaL,l) for l,deltaL in enumerate(self.phaseShifts)]
        return sum(elementsOfSum)/(self.k**2)












        

        
        


