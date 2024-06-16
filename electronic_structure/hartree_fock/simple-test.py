from Hartree_Fock import *
from representation import *
from scipy.optimize import fmin
import matplotlib.pyplot as plt
from molecular_int.GTO1s_matrix_elements import *
from molecular_int.MolecularIntegrals import *
#Find this basis, it models lithium
nullVec = np.array([0,0,0])
alphaS = [642.4180000,96.51640000,22.01740000,6.176450000,1.935110000,0.6395770000,0.5402050000,0.1022550000,0.2856450000]
typeS = [nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec]
baisS = [nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec]
alphaP = [0.5402050000,0.1022550000,0.2856450000]
basisP = [nullVec,nullVec,nullVec]
typePx = [np.array([1,0,0]),np.array([1,0,0]),np.array([1,0,0])]
typePy = [np.array([0,1,0]),np.array([0,1,0]),np.array([0,1,0])]
typePz = [np.array([0,0,1]),np.array([0,0,1]),np.array([0,0,1])]

alphaHe = [0.298073,1.242567,5.782948,38.474970]
basisHe = [nullVec,nullVec,nullVec,nullVec]
alphaH = [13.00773, 1.962079, 0.444529, 0.1219492]
basisH = [nullVec,nullVec,nullVec,nullVec]

alphas = alphaH+alphaHe+  alphaS +  alphaP + alphaP + alphaP

basisPos =basisH+ basisH+  baisS +  basisP + basisP + basisP
Zs = [3]
nucPos = [np.array([0,0,0])]
type = basisH+ basisH+typeS + typePx + typePy + typePz

rep = RepGTO(Zs,alphas,nucPos,basisPos,type)
    
t1 = [0,0,0,0,0,0,0,0]+[1 for i in range(1,19)]

t2 = [1,1,1,1,1,1,1,1]+[0 for i in range(1,19)]

ups = [t2]
downs = [t2,t1]
EGuess = -7.5
maxError = 1E-5

E,state = iterateHF(rep.normaliseList(ups),rep.normaliseList(downs),rep,EGuess,maxError,lambda s: takeGroundEigStates(s,3))
print(E)
