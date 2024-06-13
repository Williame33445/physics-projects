from Hartree_Fock import *
from representation import *

#test case from variational principle part


alphas = [0.298073,1.242567,5.782948,38.474970]


basisPos = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])]
Zs = [2]
nucPos = [np.array([0,0,0])]

rep = Rep1sGTO(Zs,alphas,nucPos,basisPos)



ups = [[1,1,1,1]]
downs = [[1,1,1,1]]
EGuess = 0
maxError = 0.1

print(iterateHF(rep.normaliseList(ups),rep.normaliseList(downs),rep,EGuess,maxError,lambda s: takeGroundEigStates(s,2)))



"""
Stuff to do:
- excited states
- hydrogen
- other cases 
- more general bases
- linear combinations of bases structure?
- make more efficient
- comments
"""