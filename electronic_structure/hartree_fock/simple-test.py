from Hartree_Fock import *
from representation import *
from scipy.optimize import fmin
import matplotlib.pyplot as plt
from molecular_int.GTO1s_matrix_elements import *
from molecular_int.MolecularIntegrals import *


def findGround(params):
    #params = [bound length,theta]
    bondLength = params[0]
    theta = params[1]
    pos = np.array([bondLength*np.cos(theta),bondLength*np.sin(theta),0])
    pos1 = np.array([bondLength,0,0])
    # 6-21G set for oxygen
    nullVec = np.array([0,0,0])
    alphaS = [0.5472270000E+04,0.8178060000E+03,0.1864460000E+03,0.5302300000E+02,0.1718000000E+02,0.5911960000E+01,0.7402940000E+01,0.1576200000E+01,0.3736840000E+00]
    typeS = [nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec]
    baisS = [nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec]

    alphaP = [0.7402940000E+01,0.1576200000E+01,0.3736840000E+00]
    basisP = [nullVec,nullVec,nullVec]
    typePx = [np.array([1,0,0]),np.array([1,0,0]),np.array([1,0,0])]
    typePy = [np.array([0,1,0]),np.array([0,1,0]),np.array([0,1,0])]
    typePz = [np.array([0,0,1]),np.array([0,0,1]),np.array([0,0,1])]


    #hydrogen
    alphasH = [13.00773,1.962079,0.444529,0.1219492]
    typeH = [nullVec,nullVec,nullVec,nullVec]
    basis1 = [pos1,pos1,pos1,pos1] #need to make naming clear
    basis2 = [pos,pos,pos,pos]



    alphas =  alphaS +  alphaP + alphaP + alphaP + alphasH + alphasH

    basisPos = baisS +  basisP + basisP + basisP + basis1 + basis2
    Zs = [8,1,1]
    nucPos = [np.array([0,0,0]),np.array([bondLength,0,0]),pos]
    type = typeS + typePx + typePy + typePz + typeH + typeH

    rep = RepGTO(Zs,alphas,nucPos,basisPos,type)
        
    t1 =[1 for i in range(1,19)] + [0 for i in range(8)]
    t2 =[0 for i in range(1,19)]  + [1,1,1,1,1,1,1,1]


    ups = [t1,t1,t1,t1,t2]
    downs = [t1,t1,t1,t1,t2]
    EGuess = 0
    maxError = 1E-4

    E,state = iterateHF(rep.normaliseList(ups),rep.normaliseList(downs),rep,EGuess,maxError,lambda s: takeGroundEigStates(s,10))
    return E + 16/bondLength + 1/(bondLength*np.sqrt(2-2*np.cos(theta)))

print(fmin(findGround,np.array([1.795,1.824])))