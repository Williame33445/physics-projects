from Hartree_Fock import *
from representation import *
from scipy.optimize import fmin
import matplotlib.pyplot as plt
from molecular_int.GTO1s_matrix_elements import *
from molecular_int.MolecularIntegrals import *
#6-21G
# nullVec = np.array([0,0,0])
# alphaS = [642.4180000,96.51640000,22.01740000,6.176450000,1.935110000,0.6395770000,0.5402050000,0.1022550000,0.2856450000]
# typeS = [nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec]
# baisS = [nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec]
# alphaP = [0.5402050000,0.1022550000,0.2856450000]
# basisP = [nullVec,nullVec,nullVec]
# typePx = [np.array([1,0,0]),np.array([1,0,0]),np.array([1,0,0])]
# typePy = [np.array([0,1,0]),np.array([0,1,0]),np.array([0,1,0])]
# typePz = [np.array([0,0,1]),np.array([0,0,1]),np.array([0,0,1])]


# # #basis is 15

# alphas = alphaS + alphaP + alphaP + alphaP

# basisPos = baisS + basisP + basisP + basisP
# Zs = [3]
# nucPos = [np.array([0,0,0])]
# type = typeS + typePx + typePy + typePz

# rep = Rep1s2p(Zs,alphas,nucPos,basisPos,type)
    


# ups = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]
# downs = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]
# EGuess = -2
# maxError = 0.01

# print(iterateHF(rep.normaliseList(ups),rep.normaliseList(downs),rep,EGuess,maxError,lambda s: takeGroundEigStates(s,3)))


# #5-21G
alphas = [0.298073,1.242567,5.782948,38.474970]
#alphas =

basisPos = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])]
Zs = [2]
nucPos = [np.array([0,0,0])]

rep1 = Rep1s2p(Zs,alphas,nucPos,basisPos,basisPos)
rep2 = Rep1sGTO(Zs,alphas,nucPos,basisPos)

print(rep1.twoElecInts-rep2.twoElecInts)

ups = [[1,1,1]]
downs = [[1,1,1]]
EGuess = -2
maxError = 0.01

print(iterateHF(rep1.normaliseList(ups),rep1.normaliseList(downs),rep1,EGuess,maxError,lambda s: takeGroundEigStates(s,2)))
print(iterateHF(rep2.normaliseList(ups),rep2.normaliseList(downs),rep2,EGuess,maxError,lambda s: takeGroundEigStates(s,2)))

# R_a = np.array([0,0,0])
# R_b = np.array([0,0,0])
# R_c = np.array([0,0,0])
# R_d = np.array([0,0,0])
# alpha_a = 0.298073
# alpha_b = 1.242567
# alpha_c = 5.782948
# alpha_d = 38.474970

# print(twoElecInt2s(R_a,alpha_a,R_b,alpha_b,R_c,alpha_c,R_d,alpha_d)- electron_repulsion(alpha_a,R_a,R_a,alpha_b,R_b,R_b,alpha_c,R_c,R_c,alpha_d,R_d,R_d))
# print(twoElecInt2s(R_a,alpha_a,R_b,alpha_b,R_c,alpha_c,R_d,alpha_d))
# print(electron_repulsion(alpha_a,R_a,R_a,alpha_b,R_b,R_b,alpha_c,R_c,R_c,alpha_d,R_d,R_d))