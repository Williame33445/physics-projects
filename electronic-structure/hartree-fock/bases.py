from abc import abstractmethod, ABC
import numpy as np 
from scipy import special as sp


#this will not work for linear combinations! May be able to double define class?

class Basis(ABC):
    def __init__(self,Zs,alphas,nuclearPositions,basisPositions):
        self.Zs = Zs
        self.alphas = alphas
        self.basisPositions = basisPositions
        self.nuclearPositions = nuclearPositions #is this the nicest way to do this?
        self.basisNumber = len(alphas)

        #find base matricies and 2 electron integrals
        self.S = self.findS()
        self.h = self.findh()
        self.twoElecInts = self.findTwoElecInts()

    def findS(self):
        S = np.empty([self.basisNumber,self.basisNumber])
        for i in range(self.basisNumber):
            for j in range(i):
                integralij = self.overLapInt(i,j)
                S[i,j] = integralij
                S[j,i] = integralij
        return S
    
    def findh(self):
        h = np.empty([self.basisNumber,self.basisNumber])
        for i in range(self.basisNumber):
            for j in range(i):
                kineticIntegral = 0.5*self.kineticInt(i,j) 
                nuclearIntegral = sum([self.nucInt(i,j,R,self.Zs[n]) for n,R in enumerate(self.nuclearPositions)])
                integral = kineticIntegral + nuclearIntegral
                h[i,j] = integral
                h[j,i] = integral
        return h
    
    def findTwoElecInts(self):
        twoElecInts = np.empty([self.basisNumber,self.basisNumber,self.basisNumber,self.basisNumber])
        for i in range(self.basisNumber):
            for k in range(i):
                for j in range(self.basisNumber):
                    for l in range(j):
                        twoElecInts[i,j,k,l] = self.twoElecInt(i,j,k,l) #this can be done with 1/2 the efficiency
        return twoElecInts

    @abstractmethod
    def overLapInt(self,p,q): #p and q are ints and specify the basis abstractly
        pass

    @abstractmethod
    def kineticInt(self,p,q):
        pass

    @abstractmethod
    def nucInt(self,p,q,R_c,Z):
        pass

    @abstractmethod
    def twoElecInt(self,p,r,q,s):
        pass



class Basis1sGTO(Basis):
    def __init__(self,Zs,alphas,nuclearPositions,basisPositions):
        Basis.__init__(self,Zs,alphas,nuclearPositions,basisPositions)


    def overLapInt(self,a,b):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        R_a = self.basisPositions[a]
        R_b = self.basisPositions[b]

        term1 = (np.pi/(alpha_a + alpha_b))**1.5
        term2 = np.exp(-(alpha_a*alpha_b)*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b))

        return term1*term2

    def kineticInt(self,a,b):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        R_a = self.basisPositions[a]
        R_b = self.basisPositions[b]

        term1 = alpha_a*alpha_b/(alpha_a + alpha_b)
        term2 = 6 - 4*alpha_a*alpha_b*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b)
        term3 = (np.pi/(alpha_a + alpha_b))**1.5
        term4 = np.exp(-(alpha_a*alpha_b)*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b))

        return term1*term2*term3*term4
    
    def nucInt(self,a,b,R_C,Z):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        R_a = self.basisPositions[a]
        R_b = self.basisPositions[b]
        R_P = (alpha_a*R_a + alpha_b*R_b)/(alpha_a + alpha_b)

        term1 = -2*np.pi*Z/(alpha_a + alpha_b)
        term2 = np.exp(-(alpha_a*alpha_b)*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b))

        t = ((alpha_a + alpha_b)*(np.linalg.norm(R_P - R_C)**2))**0.5
        
        if t == 0:
            term3 = 1
        else:
            term3 = np.sqrt(np.pi)*sp.erf(t)/(2*t)

        return term1*term2*term3

    def twoElecInt(self,a,b,c,d):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        alpha_c = self.alphas[c]
        alpha_d = self.alphas[d]
        R_a = self.basisPositions[a]
        R_b = self.basisPositions[b]
        R_c = self.basisPositions[c]
        R_d = self.basisPositions[d]

        R_P = (alpha_a*R_a + alpha_c*R_c)/(alpha_a + alpha_c)
        R_Q = (alpha_b*R_b + alpha_d*R_d)/(alpha_b + alpha_d)

        term1 = 2*(np.pi**2.5)/((alpha_a + alpha_c)*(alpha_b + alpha_d)*((alpha_a + alpha_b + alpha_c + alpha_d)**0.5))
        term2 = np.exp(-alpha_a*alpha_c*(np.linalg.norm(R_a - R_c)**2)/(alpha_a + alpha_c) - alpha_b*alpha_d*(np.linalg.norm(R_b - R_d)**2)/(alpha_b + alpha_d))
        
        
        t = np.linalg.norm(R_P - R_Q)*((alpha_a + alpha_c)*(alpha_b + alpha_d)/(alpha_a + alpha_b + alpha_c + alpha_d))**0.5
        #print(t)
        
        if t == 0:
            term3 = 1
        else:
            term3 = np.sqrt(np.pi)*sp.erf(t)/(2*t)


        # if term1*term2*term3 > 100:
        #     return 0
        return term1*term2*term3