from abc import abstractmethod, ABC
import numpy as np


class HF(ABC):
    def __init__(self,Zs,N,basis):
        self.Zs = Zs
        self.N = N #number of electrons
        self.basis = basis

    def iterateHF(self,F):
        #find ground states
        eigenvalues,eigenstates =  np.linalg.eig(np.matmul(np.linalg.inv(self.basis.S),F))
        unsortedStates = [{"epsilon": eigenvalue, "state": self.normalise(eigenstates[i])} 
                     for i, eigenvalue in enumerate(eigenvalues)]
        sortedStates = sorted(unsortedStates,key=lambda v: v["epsilon"])
        groundStates = sortedStates[-self.N:]

        #should this just be returned or should ground state energy be found



        
    def normalise(self,C):
        innerProduct = 0
        for i,c_i in enumerate(C):
            for j,c_j in enumerate(C):
                innerProduct += c_i*c_j*self.basis.S[i][j]
        return C/(innerProduct**0.5)  

    @abstractmethod
    def F(self):
        pass

class RHF(HF):
    def __init__(self,Z,basis,CGuess):
        HF.__init__(self,Z,basis)
        self.C = CGuess

    def F(self):
        pass

class UHF(HF):
    def __init__(self,Z,basis):
        HF.__init__(self,Z,basis)