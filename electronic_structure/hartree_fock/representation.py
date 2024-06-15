from abc import abstractmethod, ABC
import numpy as np 
import os
import sys
from functools import partial
from scipy import misc


from molecular_int.GTO1s_matrix_elements import *
from Hartree_Fock import *
from molecular_int.MolecularIntegrals import *





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

        return overlap1s(R_a,alpha_a,R_b,alpha_b)

    def kineticInt(self,a,b):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        R_a = self.basisPositions[a]
        R_b = self.basisPositions[b]

        return kineticInt1s(R_a,alpha_a,R_b,alpha_b)
    
    def nucInt(self,a,b,R_C,Z):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        R_a = self.basisPositions[a]
        R_b = self.basisPositions[b]

        return nucInt1s(R_a,alpha_a,R_b,alpha_b,R_C,Z)

    def twoElecInt(self,a,b,c,d):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        alpha_c = self.alphas[c]
        alpha_d = self.alphas[d]
        R_a = self.basisPositions[a]
        R_b = self.basisPositions[b]
        R_c = self.basisPositions[c]
        R_d = self.basisPositions[d]

        return twoElecInt2s(R_a,alpha_a,R_b,alpha_b,R_c,alpha_c,R_d,alpha_d)


class Rep1s2p(Representation):
    def __init__(self,Zs,alphas,nuclearPositions,basisPositions,type):
        self.alphas = alphas
        self.basisPositions = basisPositions #all R's need to be arrays
        self.type = type #[[i,j,k],...] points in direction that the integral is performed, 0 if s

        Representation.__init__(self,Zs,nuclearPositions,len(alphas))

    def overLapInt(self,a,b):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        R_a = self.basisPositions[a]
        R_b = self.basisPositions[b]
        type_a = self.type[a] #need to be tuples, also better names?
        type_b = self.type[b]     
        
        return overlap(alpha_a,type_a,R_a,alpha_b,type_b,R_b)
    
    def kineticInt(self,a,b):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        R_a = self.basisPositions[a]
        R_b = self.basisPositions[b]
        type_a = self.type[a]
        type_b = self.type[b]

        

        return 2*kinetic(alpha_a,type_a,R_a,alpha_b,type_b,R_b)
    
    def nucInt(self,a,b,R_C,Z):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        R_a = self.basisPositions[a]
        R_b = self.basisPositions[b]
        type_a = self.type[a]
        type_b = self.type[b]

        return -Z*nuclear_attraction(alpha_a,type_a,R_a,alpha_b,type_b,R_b,R_C)
    def twoElecInt(self,a,b,c,d):
        alpha_a = self.alphas[a]
        alpha_b = self.alphas[b]
        alpha_c = self.alphas[c]
        alpha_d = self.alphas[d]
        R_a = self.basisPositions[a]
        R_b = self.basisPositions[b]
        R_c = self.basisPositions[c]
        R_d = self.basisPositions[d]
        type_a = self.type[a]
        type_b = self.type[b]
        type_c = self.type[c]
        type_d = self.type[d]


        return electron_repulsion(alpha_a,type_a,R_a,alpha_b,type_b,R_b,alpha_c,type_c,R_c,alpha_d,type_d,R_d)