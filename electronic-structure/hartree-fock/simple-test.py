import sys
import os

from Hartree_Fock import *
from representation import *
from eigenstates import *

#test case from variational principle part


alphas = [0.298073,1.242567,5.782948,38.474970]

#[13.62670000,1.999350000,0.3829930000]

basisPos = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])]
Zs = [2]
nucPos = [np.array([0,0,0])]

basis = Rep1sGTO(Zs,alphas,nucPos,basisPos)


#N,basis,eigenstatesGuess,EGuess,maxError

eigenstatesGuess = Eigenstates([[1,1,1,1]],[[1,1,1,1]],basis)
EGuess = 0
maxError = 0.01

#print(np.linalg.inv(basis.S))
HF(Zs,2,basis,eigenstatesGuess,EGuess,maxError)


#should look at hydrogen, put this case into a test, make more efficient and develop the class more
# excited states?