from abc import abstractmethod, ABC
import numpy as np 

class RHF(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def overLapInt(self,p,q): #p and q are ints and specify the basis abstractly
        pass

    @abstractmethod
    def kineticInt(self,p,q):
        pass

    @abstractmethod
    def nucInt(self,p,q):
        pass

    @abstractmethod
    def twoElecInt(self,p,r,q,s):
        pass



class RHF1sGTOBasis(RHF):
    def __init__(self,alphas,nuclearPositions):
        #all abstract class parts are defined first
        RHF.__init__(self)
        self.alphas = alphas
        self.nuclearPositions = nuclearPositions 


    def overLapInt(self,p,q):
        alpha_p = self.alphas[p]
        alpha_q = self.alphas[q]
        R_p = self.nuclearPositions[p]
        R_q = self.nuclearPositions[q]

        term1 = (np.pi/(alpha_p + alpha_q))**1.5
        term2 = np.exp(-(alpha_p*alpha_q)*(np.linalg.norm(R_p - R_q)**2)/(alpha_p + alpha_q))

        return term1*term2

    def kineticInt(self,p,q):
        alpha_p = self.alphas[p]
        alpha_q = self.alphas[q]
        R_p = self.nuclearPositions[p]
        R_q = self.nuclearPositions[q]

        term1 = alpha_p*alpha_q/(alpha_p + alpha_q)
        term2 = 6 - 4*alpha_p*alpha_q*(np.linalg.norm(R_p - R_q)**2)/(alpha_p + alpha_q)
        term3 = (np.pi/(alpha_p + alpha_q))**1.5
        term4 = np.exp(-(alpha_p*alpha_q)*(np.linalg.norm(R_p - R_q)**2)/(alpha_p + alpha_q))

        return term1*term2*term3*term4
    
    def nucInt(self,p,q):
        pass

    def twoElecInt(self,p,r,q,s):
        pass
