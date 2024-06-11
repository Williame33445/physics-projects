from abc import abstractmethod, ABC
import numpy as np 
from scipy import special as sp


class Basis(ABC):
    def __init__(self,Z,alphas,nuclearPositions):
        self.Z = Z
        self.alphas = alphas
        self.nuclearPositions = nuclearPositions 

    @abstractmethod
    def overLapInt(self,p,q): #p and q are ints and specify the basis abstractly
        pass

    @abstractmethod
    def kineticInt(self,p,q):
        pass

    @abstractmethod
    def nucInt(self,p,q,R_c):
        pass

    @abstractmethod
    def twoElecInt(self,p,r,q,s):
        pass



class Basis1sGTO(Basis):
    def __init__(self,Z,alphas,nuclearPositions):
        Basis.__init__(self,Z,alphas,nuclearPositions)


    def overLapInt(self,a,b):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        R_a = self.nuclearPositions[a]
        R_b = self.nuclearPositions[b]

        term1 = (np.pi/(alpha_a + alpha_b))**1.5
        term2 = np.exp(-(alpha_a*alpha_b)*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b))

        return term1*term2

    def kineticInt(self,a,b):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        R_a = self.nuclearPositions[a]
        R_b = self.nuclearPositions[b]

        term1 = alpha_a*alpha_b/(alpha_a + alpha_b)
        term2 = 6 - 4*alpha_a*alpha_b*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b)
        term3 = (np.pi/(alpha_a + alpha_b))**1.5
        term4 = np.exp(-(alpha_a*alpha_b)*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b))

        return term1*term2*term3*term4
    
    def nucInt(self,a,b,R_C):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        R_a = self.nuclearPositions[a]
        R_b = self.nuclearPositions[b]
        R_P = (alpha_a*R_a + alpha_b*R_b)/(alpha_a + alpha_b)

        term1 = -2*np.pi*self.Z/(alpha_a + alpha_b)
        term2 = np.exp(-(alpha_a*alpha_b)*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b))

        t = ((alpha_a + alpha_b)*(np.linalg.norm(R_P - R_C)**2))**0.5
        term3 = np.sqrt(np.pi)*np.erf(t)/(2*t)

        return term1*term2*term3

    def twoElecInt(self,a,b,c,d):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        alpha_c = self.alphas[c]
        alpha_d = self.alphas[d]
        R_a = self.nuclearPositions[a]
        R_b = self.nuclearPositions[b]
        R_c = self.nuclearPositions[c]
        R_d = self.nuclearPositions[d]

        R_P = (alpha_a*R_a + alpha_c*R_c)/(alpha_a + alpha_c)
        R_Q = (alpha_b*R_b + alpha_d*R_d)/(alpha_b + alpha_d)

        term1 = 2*(np.pi**2.5)/((alpha_a + alpha_c)*(alpha_b + alpha_d)*((alpha_a + alpha_b + alpha_c + alpha_d)**0.5))
        term2 = np.exp(-alpha_a*alpha_c*(np.linalg.norm(R_a - R_c)**2)/(alpha_a + alpha_c) - alpha_b*alpha_d*(np.linalg.norm(R_b - R_d)**2)/(alpha_b + alpha_d))
        
        t = np.linalg.norm(R_P - R_Q)*((alpha_a + alpha_c)*(alpha_b + alpha_d)/(alpha_a + alpha_b + alpha_c + alpha_d))**0.5
        term3 = np.sqrt(np.pi)*np.erf(t)/(2*t)

        return term1*term2*term3