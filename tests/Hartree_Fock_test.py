import unittest
import os
import sys
import numpy as np

sys.path.append(os.path.abspath("."))

from electronic_structure.hartree_fock.Hartree_Fock import *
from electronic_structure.hartree_fock.representation import *


class Test_Hartree_Fock(unittest.TestCase):

    def testHeliumGround(self):
        alphas = [0.298073,1.242567,5.782948,38.474970]
        basisPos = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])]
        Zs = [2]
        nucPos = [np.array([0,0,0])]

        rep = Rep1sGTO(Zs,alphas,nucPos,basisPos)
        ups = [[1,1,1,1]]
        downs = [[1,1,1,1]]
        EGuess = 0
        maxError = 0.1

        result = iterateHF(rep.normaliseList(ups),rep.normaliseList(downs),rep,EGuess,maxError,lambda s: takeGroundEigStates(s,2))
        expected = -2.81

        difference = np.round(result - expected, decimals=2)
            
        self.assertEqual(0,difference)
        