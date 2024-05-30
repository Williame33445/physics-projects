import unittest
import os
import sys
import numpy as np
from spherically_symetric_quantum_scattering import ScatteringSystem

class TestNumerov(unittest.TestCase):

    def testQHO(self):
        #system parameters
        V = lambda r : r**2
        E = 3
        l = 0
        m = 1

        #simulation parameters
        h = 0.1

        #intial conditions
        r_0 = 0
        r_end = 4
        u_0 = 0
        u_1 = h**(l+1) #this is a guess from solving a well behaved system at the origin

        
        
        scatteringSys = ScatteringSystem(E,l,m,V,r_0,h,u_0,u_1,r_end,simpleUnits=True)

        #finds expected predictions
        rArray = np.array(scatteringSys.rList)
        uExpected = np.exp((h**2)/2)*rArray*np.exp(-(rArray**2)/2)
   

        for i,u in enumerate(scatteringSys.uList):
            #as h^4 is the order than this method is correct to at least 2 dp
            #allows for some fluctuations due to rounding

            difference = np.round(u - uExpected[i], decimals=2)
            
            self.assertEqual(0,difference)


