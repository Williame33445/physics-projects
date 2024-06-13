import sys
import os

from Hartree_Fock import *
from representation import *

#test case from variational principle part


alphas = [0.298073,1.242567,5.782948,38.474970]


basisPos = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])]
Zs = [2]
nucPos = [np.array([0,0,0])]

rep = Rep1sGTO(Zs,alphas,nucPos,basisPos)


#N,basis,eigenstatesGuess,EGuess,maxError 

ups = [[1,1,1,1]]
downs = [[1,1,1,1]]
EGuess = 0
maxError = 0.1

#print(np.linalg.inv(basis.S))
print(iterateHF(rep.normaliseList(ups),rep.normaliseList(downs),rep,2,EGuess,maxError))


#should look at hydrogen, put this case into a test, make more efficient and develop the class more
# excited states?