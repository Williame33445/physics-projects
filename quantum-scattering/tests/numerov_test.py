import unittest
import os
import sys
import numpy as np

sys.path.append(os.path.abspath("./quantum-scattering"))

from numerov_method import *



#setup for QHO test
def F(l,r,E,V):
    if r == 0:
        #allowed to do this as it only occurs at u_0 = 0, so it has no effect (make this clearer)
        return 0
    return V(r) + l*(l+1)/(r**2) - E

#setup for hydrogen test
def f(E):
    return lambda r: -2*(E + 1/r)

def uAsym(r):
    return r*np.exp(-r) 

def getRMin(E):
    rs,us = runNumerov(100,-0.1,uAsym(100),uAsym(100-0.1),f(E),999)
    return us[-1]




class TestNumerov(unittest.TestCase):

    def testQHO(self):
        """
        This is a simple implementation of the numerov test outlined in chapter 3 in the computational physics book
        """
        #simulation setup
        V = lambda r : r**2
        E = 3
        l = 0
        FQHO1stEx = lambda r: F(l,r,E,V)

        #simulation parameters
        h = 0.1
        N = 30
        u_0 = 0
        u_1 = h**(l+1)
        
        #run simulation
        rList,uNumericalList = runNumerov(0,h,u_0,u_1,FQHO1stEx,N)

        #finds expected predictions
        rArray = np.array(rList)
        uExpected = np.exp((h**2)/2)*rArray*np.exp(-(rArray**2)/2)        

        for i,u in enumerate(uNumericalList):
            #as h^4 is the order than this method is correct to at least 2 dp
            #allows for some fluctuations due to rounding

            difference = np.round(u - uExpected[i], decimals=2)
            
            self.assertEqual(0,difference)
    
    def testHydrogen(self):
        """
        This is a simple method for finding the hydrogen atoms ground state.
        """
        

        ERange = -np.arange(0.01,2,0.01)
        EGround = min(ERange,key=getRMin)
        EExpected = -0.5

        difference = np.round(EGround - EExpected, decimals=1)
        self.assertEqual(0,difference)




if __name__=="__main__":
    unittest.main()