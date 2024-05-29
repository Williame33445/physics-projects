import unittest
import os
import sys
import numpy as np

sys.path.append(os.path.abspath("."))

from useful_code.numerov_method import *
from hydrogen_scattering_by_krypton_atoms import * 

class TestNumerov(unittest.TestCase):

    def testQHO(self):
        #simulation setup
        V = lambda r : r**2
        E = 3
        l = 0
        FQHO1stEx = lambda r: F(l,r,E,1,V,units=True)

        #simulation parameters
        h = 0.1
        N = 30
        u_0 = 0
        u_1 = h**(l+1) #this is a guess from solving a well behaved system at the origin
        
        #run simulation
        rList,uNumericalList = runNumerov(u_0,u_1,h,FQHO1stEx,N)

        #finds expected predictions
        rArray = np.array(rList)
        uExpected = np.exp((h**2)/2)*rArray*np.exp(-(rArray**2)/2)        

        for i,u in enumerate(uNumericalList):
            #as h^4 is the order than this method is correct to at least 2 dp
            #allows for some fluctuations due to rounding

            difference = np.round(u - uExpected[i], decimals=2)
            
            self.assertEqual(0,difference)


if __name__=="__main__":
    unittest.main()


