import sys
import os
import numpy as np
from scipy import special

sys.path.append(os.path.abspath("."))
from useful_code.numerov_method import *

#units, defines 2m/hbar^2 units
units = 6.12

class ScatteringSystem:
    """
    Describes an instance of scattering of two molecules that interact under the Lennard-Jones potential.
    Uses partial wave analysis and the Numerov method. 

    Inputs:
    E - energy of the system
    V - spherically symmetric potential of the system
    r_0 - first initial condition r value
    h - distance between r_0 and the r value corresponding to u_1
    u_0 - first initial condition u value
    u_1 - second initial condition u value
    r_max - boundary radius for the V != 0 region
    lCutoff - value of l that the sums are performed to (usually 9)
    """
    def __init__(self,E,V,r_0,h,u_0,u_1,r_max,lCutoff):
        #system parameters
        self.E = E
        self.k = np.sqrt(units*E)
        self.wavelength = 2*np.pi/self.k
        self.V = V
        self.lCutOff = lCutoff
        self.h = h
        self.r_0 = r_0

        #finds max number of points separated by h required to describe (r < b1 + wavelength) region
        self.N = int(np.ceil((r_max + self.wavelength - r_0)/h)) 

        #finds b1 index and the r corresponding to that index, used as boundary conditions for finding phases 
        self.r_b1Index = int(np.ceil((r_max-r_0)/h))
        self.r_b1 = h*self.r_b1Index + r_0 

        #finds final r value calculated, used as boundary condition for finding phases
        self.r_b2 = r_0 + self.N*h

        #runs simulations for l values
        self.uLists = []
        for l in range(self.lCutOff):
            F = lambda r: self.FL(r,l)
            self.rList, ulList = runNumerov(r_0,h,u_0,u_1,F,self.N)
            self.uLists.append(ulList)

        #finds phase shifts and cross sections
        self.phaseShifts = self.findPhaseShifts()
        self.totalCrossSection = self.findTotalCrossSection()

    def FL(self,r,l):
        """
        The F_l(r) part of the Numerov method.
        """
        return units*self.V(r) + l*(l+1)/(r**2) - units*self.E
    
    def findlPhaseShift(self,l,u_max,u_end):
        """
        Given boundary conditions finds l phase shift from eq (5).
        """
        K = (self.r_b1*u_end)/(self.r_b2*u_max)
        numerator = K*special.spherical_jn(l,self.k*self.r_b1) - special.spherical_jn(l,self.k*self.r_b2)
        denominator = K*special.spherical_yn(l,self.k*self.r_b1) - special.spherical_yn(l,self.k*self.r_b2)
        return np.arctan(numerator/denominator)


    def findPhaseShifts(self):
        """
        Finds all phase shifts from uLists.
        """
        phaseShifts = []
        for l in range(self.lCutOff):
            u_max = self.uLists[l][self.r_b1Index]
            u_end = self.uLists[l][-1] #as u_end at the size of the list
            phaseShifts.append(self.findlPhaseShift(l,u_max,u_end))
        return phaseShifts
            
    def findTotalCrossSection(self):
        """
        Uses eq (7) to find total cross section.
        """
        elementsOfSum = [(2*l + 1)*(np.sin(deltaL)**2) for l,deltaL in enumerate(self.phaseShifts)]
        return 4*np.pi*sum(elementsOfSum)/(self.k**2)

    def differentialCrossSection(self,theta):
        """
        Uses eq (6) to find differential cross section.
        """
        termOfSum = lambda deltaL,l:  np.abs((2*l+1)*np.exp(complex(0,1)*deltaL)*np.sin(deltaL)*
                                special.eval_legendre(l,np.cos(theta)))**2
        elementsOfSum = [termOfSum(deltaL,l) for l,deltaL in enumerate(self.phaseShifts)]
        return sum(elementsOfSum)/(self.k**2)
