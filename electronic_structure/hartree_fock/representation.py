from abc import abstractmethod, ABC
import numpy as np 
from scipy import special as sp

def F_0(t):
    if t == 0:
        return 1
    else:
        return np.sqrt(np.pi)*sp.erf(t)/(2*t)



#can this be generalised for linear combinations two bases?
#holds and describes the representations of the matricies

class Representation(ABC):
    def __init__(self,Zs,nuclearPositions,repNumb):
        self.Zs = Zs
        self.nuclearPositions = nuclearPositions 
        self.repNumb = repNumb

        #find base matricies and 2 electron integrals
        self.S = self.findS()
        self.h = self.findh()
        self.twoElecInts = self.findTwoElecInts()

        self.invS = np.linalg.inv(self.S)

    def findS(self):
        S = np.empty([self.repNumb,self.repNumb])
        for i in range(self.repNumb):
            for j in range(i+1):
                integralij = self.overLapInt(i,j)
                S[i,j] = integralij
                S[j,i] = integralij
        return S
    
    def findh(self):
        h = np.empty([self.repNumb,self.repNumb])
        for i in range(self.repNumb):
            for j in range(i+1):
                kineticIntegral = 0.5*self.kineticInt(i,j) 
                nuclearIntegral = sum([self.nucInt(i,j,R,self.Zs[n]) for n,R in enumerate(self.nuclearPositions)])
                integral = kineticIntegral + nuclearIntegral
                h[i,j] = integral
                h[j,i] = integral
        return h
    
    def findTwoElecInts(self):
        twoElecInts = np.empty([self.repNumb,self.repNumb,self.repNumb,self.repNumb])
        for i in range(self.repNumb):
            for k in range(self.repNumb):
                for j in range(self.repNumb):
                    for l in range(self.repNumb):
                        twoElecInts[i,j,k,l] = self.twoElecInt(i,j,k,l) #this can be done with 1/2 the efficiency
        return twoElecInts
    
    def normalise(self,C):
        innerProduct = 0
        for i,c_i in enumerate(C):
            for j,c_j in enumerate(C):
                innerProduct += c_i*c_j*self.S[i][j]
        return C/(innerProduct**0.5)
    
    def normaliseList(self,lst):
        for i,C in enumerate(lst):
            lst[i] = self.normalise(C)
        return lst
    
    def normaliseDicList(self,lst):
        for i,D in enumerate(lst):
            lst[i]["state"] = self.normalise(D["state"])
        return lst
    
    def F(self,states1,states2):
        F = np.empty([self.repNumb,self.repNumb])
        for p in range(self.repNumb):
            for q in range(self.repNumb):
                F[p,q] = self.h[p,q]

                for C in states1:
                    for r,c_r in enumerate(C): #can be made more efficient
                        for s,c_s in enumerate(C):
                            F[p,q] += c_r*c_s*(self.twoElecInts[p,r,q,s] - self.twoElecInts[p,r,s,q])

                
                for C in states2:
                    for r,c_r in enumerate(C):
                        for s,c_s in enumerate(C): 
                            F[p,q] += c_r*c_s*self.twoElecInts[p,r,q,s]
        
        return F
    
    
    def findE(self,states):
        E = 0
        for S in states:
            E += S["e"]
            for i,c_i in enumerate(S["state"]):
                for j,c_j in enumerate(S["state"]):
                    E += self.h[i,j]*c_i*c_j
        
        return E/2


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



class Rep1sGTO(Representation):
    def __init__(self,Zs,alphas,nuclearPositions,basisPositions):
        self.alphas = alphas
        self.basisPositions = basisPositions
        Representation.__init__(self,Zs,nuclearPositions,len(alphas))
        
        


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
        
        term3 = F_0(t)

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
        
        term3 = F_0(t)

        return term1*term2*term3
