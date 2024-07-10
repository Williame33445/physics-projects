# performs hartree fock calculation on a cooper atom
# 3-21G basis from bse

import numpy as np
import sys
import os
import json

sys.path.append(os.path.abspath("."))

from electronic_structure.hartree_fock.Hartree_Fock import *
from electronic_structure.hartree_fock.representation import *

#sets up basis
nullVec = np.array([0,0,0])
alphaS = [0.4134302200E+04,0.6254912200E+03,0.1369555600E+03,
          0.1814960330E+03,0.3957431190E+02,0.1216246380E+02,
          0.1235111490E+02,0.4049651020E+01,0.1279225380E+01,
          0.1048299940E+01,0.1171180220E+00,
          0.4054498690E-01]
typeS = [nullVec for i in range(12)]
baisS = [nullVec for i in range(12)]

alphaP = [0.1814960330E+03,0.3957431190E+02,0.1216246380E+02,
          0.1235111490E+02,0.4049651020E+01,0.1279225380E+01,
          0.1048299940E+01,0.1171180220E+00,
          0.4054498690E-01]
basisP = [nullVec for i in range(9)]
typePx = [np.array([1,0,0]) for i in range(9)]
typePy = [np.array([0,1,0]) for i in range(9)]
typePz = [np.array([0,0,1]) for i in range(9)]

alphaD = [0.1675937600E+02,0.4178976900E+01,0.9943270400E+00]
basisD = [nullVec for i in range(3)]
typeDxx = [np.array([2,0,0]) for i in range(3)]
typeDyy = [np.array([0,2,0]) for i in range(3)]
typeDzz = [np.array([0,0,2]) for i in range(3)]
typeDxy = [np.array([1,1,0]) for i in range(3)]
typeDxz = [np.array([1,0,1]) for i in range(3)]
typeDyz = [np.array([0,1,1]) for i in range(3)]

alphas = alphaS + alphaP + alphaP + alphaP + alphaD + alphaD + alphaD + alphaD + alphaD + alphaD
basisPos= baisS + basisP + basisP + basisP + basisD + basisD + basisD + basisD + basisD + basisD
type =    typeS + typePx + typePy + typePz + typeDxx + typeDyy + typeDzz + typeDxy + typeDxz + typeDyz

#sets up simulation parameters
Zs = [29]
nucPos = [np.array([0,0,0])]
maxError = 1E-4

#sets up representation
rep = RepGTO(Zs,alphas,nucPos,basisPos,type)

print("found integrals")

#sets up intial guess  
g =[1 for i in range(57)]
ups = [g for i in range(14)]
downs = [g for i in range(15)]
EGuess = 0

#finds energy and prints
E,state = iterateHF(rep.normaliseList(ups),rep.normaliseList(downs),rep,EGuess,maxError,lambda s: takeGroundEigStates(s,8))
print(f"Energy: {E} (hartree)")

with open("electronic_structure\periodic-solids\copper.json", "w") as f:
    json.dump(E,f)
    json.dump(state,f)