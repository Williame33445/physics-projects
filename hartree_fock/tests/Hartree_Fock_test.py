import unittest
import os
import sys
import numpy as np

sys.path.append(os.path.abspath("./hartree_fock"))

from Hartree_Fock import *
from representation import *
from get_GTOs_from_BSE import *
from GTOs import *

class Test_Hartree_Fock(unittest.TestCase):

    def testHeliumGround(self):
        #sets up primitive 3ZaPa-NR-CV Helium basis set (should be better for excited states)
        GTOsHe = getGuassians("3ZaPa-NR-CV", 2, np.zeros(3), basisType="primitive")

        #simulation parameters
        ZsHe = [2]
        nucPosHe = [np.array([0, 0, 0])]
        maxErrorHe = 1E-4

        #define representation class
        repHe = RepGTO(GTOsHe, ZsHe, nucPosHe)

        #define intial guess
        ups = repHe.normaliseList([[1 for i in range(len(GTOsHe))]])
        downs = repHe.normaliseList([[1 for i in range(len(GTOsHe))]])
        EGuess = -2.8

        #calculate ground state and print energy
        E1s1s_S, HeElec1s1s_S = iterateHF(ups, downs, repHe, EGuess, maxErrorHe, lambda s: takeGroundEigStates(s,2))

        expected = -2.86158
        difference = np.round(E1s1s_S - expected, decimals=1)
            
        self.assertEqual(0,difference)

if __name__=="__main__":
    unittest.main()