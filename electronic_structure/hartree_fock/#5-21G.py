from Hartree_Fock import *
from representation import *
from scipy.optimize import fmin
import matplotlib.pyplot as plt
from molecular_int.GTO1s_matrix_elements import *
from molecular_int.MolecularIntegrals import *

#5-21G
alphas = [0.298073,1.242567,5.782948,38.474970]

basisPos = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])]
Zs = [2]
nucPos = [np.array([0,0,0])]

rep = RepGTO(Zs,alphas,nucPos,basisPos,basisPos)



ups = [[1,1,1,1]]
downs = [[1,1,1,1]]
EGuess = -2.8
maxError = 1E-14

print(iterateHF(rep.normaliseList(ups),rep.normaliseList(downs),rep,EGuess,maxError,lambda s: takeGroundEigStates(s,2)))

