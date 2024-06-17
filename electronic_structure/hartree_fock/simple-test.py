from Hartree_Fock import *
from representation import *
from scipy.optimize import fmin
import matplotlib.pyplot as plt
from molecular_int.GTO1s_matrix_elements import *
from molecular_int.MolecularIntegrals import *
alphasH2 = [13.00773,1.962079,0.444529,0.1219492,13.00773,1.962079,0.444529,0.1219492]
bondLength = 1.38810547

basisPos = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,bondLength]),np.array([0,0,bondLength]),np.array([0,0,bondLength]),np.array([0,0,bondLength])]
Zs = [1,1]
nucPos = [np.array([0,0,0]),np.array([0,0,bondLength])]
alphas = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])]

rep = RepGTO(Zs,alphasH2,nucPos,basisPos,alphas)



ups = [[1,1,1,1,1,1,1,1]]
downs = [[1,1,1,1,1,1,1,1]]
EGuess = 0  
maxError = 1E-4

E_hf,states = iterateHF(rep.normaliseList(ups),rep.normaliseList(downs),rep,EGuess,maxError,lambda s: takeGroundEigStates(s,2))
#koopmans theorm for ionisation energy
print(states)