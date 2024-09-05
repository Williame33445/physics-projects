from abc import abstractmethod, ABC
import numpy as np 

from GTO1s_matrix_elements import *
from MolecularIntegrals import *


class Representation(ABC):
    """
    Abstract class that defines the general form of the representation.

    Inputs:
    Zs - list of nuclear numbers
    nuclearPositions - list of nuclear positions
    repNumb - number of functions in the basis
    """

    def __init__(self,Zs,nuclearPositions,repNumb):
        self.Zs = Zs
        self.nuclearPositions = nuclearPositions 
        self.repNumb = repNumb

        #find base matricies and 2 electron integrals
        self.S = self.findS()
        self.h = self.findh()
        self.twoElecInts = self.findTwoElecInts()


    def findS(self):
        """
        Method that finds S matrix
        """
        S = np.empty([self.repNumb,self.repNumb])
        for i in range(self.repNumb):
            for j in range(i+1):
                integralij = self.overLapInt(i,j)
                S[i,j] = integralij
                S[j,i] = integralij
        return S
    
    def findh(self):
        """
        Method that finds h matrix.
        """
        h = np.empty([self.repNumb,self.repNumb])
        for i in range(self.repNumb):
            for j in range(i+1):
                kineticIntegral = self.kineticInt(i,j) 
                nuclearIntegral = sum([self.nucInt(i,j,R,self.Zs[n]) for n,R in enumerate(self.nuclearPositions)])
                integral = kineticIntegral + nuclearIntegral
                h[i,j] = integral
                h[j,i] = integral
        return h
    
    def findTwoElecInts(self):
        """
        Method that finds 2 electron integrals.
        """
        twoElecInts = np.empty([self.repNumb,self.repNumb,self.repNumb,self.repNumb])
        for i in range(self.repNumb):
            for k in range(i+1):
                for j in range(self.repNumb):
                    for l in range(j+1):
                        integral = self.twoElecInt(i,j,k,l)
                        twoElecInts[i,j,k,l] = integral
                        twoElecInts[k,j,i,l] = integral
                        twoElecInts[i,l,k,j] = integral
                        twoElecInts[k,l,i,j] = integral
        return twoElecInts
    
    def normalise(self,C):
        """
        Method that normalises numpy array C.
        """
        return C/(np.einsum("ij,i,j",self.S,C,C)**0.5)
    
    def normaliseList(self,lst):
        """
        Method that normalises a list of numpy arrays.
        """
        return [self.normalise(C) for C in lst]
    
    def normaliseDicList(self,lst):
        """
        Method that normalises a list of dictionaries in the form 
        defined in hartree_fock.py.
        """
        for i,D in enumerate(lst):
            lst[i]["state"] = self.normalise(D["state"])
        return lst
    
    def F(self,states1,states2):
        """
        Method that finds F^(+/-) matrix given C^(+/-) and C^(-/+) matricies.
        """
        sum1 = sum([np.einsum("prqs,r,s",self.twoElecInts, C, C) for C in states1])
        sum2 = sum([np.einsum("prsq,r,s",self.twoElecInts, C, C) for C in states1])
        sum3 = sum([np.einsum("prqs,r,s",self.twoElecInts, C, C) for C in states2])
        return self.h + sum1 - sum2 +sum3
    
    
    def findE(self,states):
        """
        Method that finds E given a list of dictionaries in the form 
        defined in hartree_fock.py
        """
        return sum([S["e"] + np.einsum("ij,i,j",self.h,S["state"],S["state"]) for S in states])/2 #use S or C consistently


    @abstractmethod
    def overLapInt(self,p,q):
        """
        Abstract method that finds the overlap integral elements.
        """
        pass

    @abstractmethod
    def kineticInt(self,p,q):
        """
        Abstract method that finds the kinetic integral elements.
        """
        pass

    @abstractmethod
    def nucInt(self,p,q,R_c,Z):
        """
        Abstract method that finds the nuclear integral elements.
        """
        pass

    @abstractmethod
    def twoElecInt(self,p,r,q,s):
        """
        Abstract method that finds the two electron integral elements.
        """
        pass



class Rep1sGTO(Representation):
    """
    Representation class for the symmetric GTO basis. 
    """
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
    """
    Representation class for the general GTO basis. 
    """
    def __init__(self,Zs,alphas,nuclearPositions,basisVecs,type):
        self.alphas = alphas
        self.basisVecs = basisVecs 
        self.type = type 
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
        #order of c and b's here is due to notational differences
        return electron_repulsion(self.alphas[a],self.type[a],self.basisVecs[a],
                                  self.alphas[c],self.type[c],self.basisVecs[c],
                                  self.alphas[b],self.type[b],self.basisVecs[b],
                                  self.alphas[d],self.type[d],self.basisVecs[d])
