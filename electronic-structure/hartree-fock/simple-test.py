import sys
import os

from Hartree_Fock import *
from bases import *

#from Basis set exchange Basis set: 5-21G


alphas = [0.298073,1.242567,5.782948,38.474970]

#[13.62670000,1.999350000,0.3829930000]

basisPos = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])]
Zs = [2]
nucPos = [np.array([0,0,0])]

basis = Basis1sGTO(Zs,alphas,nucPos,basisPos)


#N,basis,eigenstatesGuess,EGuess,maxError

eigenstatesGuess = Eigenstates([[1,1,1,1]],[[1,1,1,1]])
EGuess = -2
maxError = 1

#print(np.linalg.inv(basis.S))
HF(Zs,2,basis,eigenstatesGuess,EGuess,maxError)