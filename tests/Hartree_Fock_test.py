import unittest
import os
import sys
import numpy as np

sys.path.append(os.path.abspath("."))

from electronic_structure.hartree_fock.Hartree_Fock import *
from electronic_structure.hartree_fock.representation import *


class Test_Hartree_Fock(unittest.TestCase):

    def testHeliumGround(self):
        #define simulation parameters
        alphas = [0.298073,1.242567,5.782948,38.474970]
        basisPos = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])]
        Zs = [2]
        nucPos = [np.array([0,0,0])]
        maxError = 1E-4

        #define representation class
        rep = RepGTO(Zs,alphas,nucPos,basisPos,basisPos)

        #define intial guess
        ups = [[1,1,1,1]]
        downs = [[1,1,1,1]]
        EGuess = -2.8
        expected = -2.8551714954912644

        #run simulation and print energy
        E,states = iterateHF(rep.normaliseList(ups),rep.normaliseList(downs),rep,EGuess,maxError,lambda s: takeGroundEigStates(s,2))

        difference = np.round(E - expected, decimals=1)
            
        self.assertEqual(0,difference)
        