from abc import abstractmethod, ABC
import numpy as np 


from molecular_int.GTO1s_matrix_elements import *
from Hartree_Fock import *
from molecular_int.MolecularIntegrals import *

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


    def findS(self):
        S = np.empty([self.repNumb,self.repNumb])
        for i in range(self.repNumb):
            for j in range(self.repNumb):
                integralij = self.overLapInt(i,j)
                S[i,j] = integralij
        return S
    
    def findh(self):
        h = np.empty([self.repNumb,self.repNumb])
        for i in range(self.repNumb):
            for j in range(self.repNumb):
                kineticIntegral = self.kineticInt(i,j) 
                nuclearIntegral = sum([self.nucInt(i,j,R,self.Zs[n]) for n,R in enumerate(self.nuclearPositions)])
                integral = kineticIntegral + nuclearIntegral
                h[i,j] = integral
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
        return C/(np.einsum("ij,i,j",self.S,C,C)**0.5)
    
    def normaliseList(self,lst):
        return [self.normalise(C) for C in lst]
    
    def normaliseDicList(self,lst):
        for i,D in enumerate(lst):
            lst[i]["state"] = self.normalise(D["state"])
        return lst
    
    def F(self,states1,states2):
        sum1 = sum([np.einsum("prqs,r,s",self.twoElecInts, C, C) for C in states1])
        sum2 = sum([np.einsum("prsq,r,s",self.twoElecInts, C, C) for C in states1])
        sum3 = sum([np.einsum("prqs,r,s",self.twoElecInts, C, C) for C in states2])
        return self.h + sum1 - sum2 +sum3
    
    
    def findE(self,states):
        return sum([S["e"] + np.einsum("ij,i,j",self.h,S["state"],S["state"]) for S in states])/2 #use S or C consistently


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
    def __init__(self,Zs,alphas,nuclearPositions,basisVecs):
        self.alphas = alphas
        self.basisVecs = basisVecs
        Representation.__init__(self,Zs,nuclearPositions,len(alphas))
        
    def overLapInt(self,a,b):
        return overlap1s(self.basisVecs[a],self.alphas[a],self.basisVecs[b],self.alphas[b])

    def kineticInt(self,a,b):
        return kineticInt1s(self.basisVecs[a],self.alphas[a],self.alphas[b],self.alphas[b])
    
    def nucInt(self,a,b,R_C,Z):
        return nucInt1s(self.basisVecs[a],self.alphas[a],self.alphas[b],self.alphas[b],R_C,Z)

    def twoElecInt(self,a,b,c,d):
        return twoElecInt2s(self.basisVecs[a],self.alphas[a],self.alphas[b],self.alphas[b],
                            self.basisVecs[c],self.alphas[c],self.alphas[d],self.alphas[d])


class RepGTO(Representation):
    def __init__(self,Zs,alphas,nuclearPositions,basisVecs,type):
        self.alphas = alphas
        self.basisVecs = basisVecs #all R's need to be arrays
        self.type = type #[[i,j,k],...] points in direction that the integral is performed, 0 if s
        Representation.__init__(self,Zs,nuclearPositions,len(alphas))

    def overLapInt(self,a,b):  
        return overlap(self.alphas[a],self.type[a],self.basisVecs[a],
                       self.alphas[b],self.type[b],self.basisVecs[b])
    
    def kineticInt(self,a,b):
        return kinetic(self.alphas[a],self.type[a],self.basisVecs[a],
                       self.alphas[b],self.type[b],self.basisVecs[b])
    
    def nucInt(self,a,b,R_C,Z):
        return -Z*nuclear_attraction(self.alphas[a],self.type[a],self.basisVecs[a],
                                     self.alphas[b],self.type[b],self.basisVecs[b],R_C)
    def twoElecInt(self,a,b,c,d):
        return electron_repulsion(self.alphas[a],self.type[a],self.basisVecs[a],
                                  self.alphas[c],self.type[c],self.basisVecs[c],
                                  self.alphas[b],self.type[b],self.basisVecs[b],
                                  self.alphas[d],self.type[d],self.basisVecs[d])
