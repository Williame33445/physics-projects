from Hartree_Fock import *
from representation import *
from scipy.optimize import fmin
import matplotlib.pyplot as plt
#6-21G
nullVec = np.array([0,0,0])
alphaS = [642.4180000,96.51640000,22.01740000,6.176450000,1.935110000,0.6395770000,0.5402050000,0.1022550000,0.2856450000]
typeS = [nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec]
baisS = [nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec]
alphaP = [0.5402050000,0.1022550000,0.2856450000]
basisP = [nullVec,nullVec,nullVec]
typePx = [np.array([1,0,0]),np.array([1,0,0]),np.array([1,0,0])]
typePy = [np.array([0,1,0]),np.array([0,1,0]),np.array([0,1,0])]
typePz = [np.array([0,0,1]),np.array([0,0,1]),np.array([0,0,1])]


#basis is 15

alphas = alphaS + alphaP + alphaP + alphaP

basisPos = baisS + basisP + basisP + basisP
Zs = [3]
nucPos = [np.array([0,0,0])]
type = typeS + typePx + typePy + typePz

rep = Rep1s2p(Zs,alphas,nucPos,basisPos,type)



ups = [[1,1,1,1,1,1,1,1,1,1,0,0,0,0,0]]
downs = [[1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]
EGuess = -2
maxError = 0.01

print(iterateHF(rep.normaliseList(ups),rep.normaliseList(downs),rep,EGuess,maxError,lambda s: takeGroundEigStates(s,3)))

