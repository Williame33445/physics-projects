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
    def __init__(self,GTOs,Zs,nuclearPositions):
        self.GTOs = GTOs
        Representation.__init__(self,Zs,nuclearPositions,len(GTOs))
        
    def overLapInt(self,a,b):
        A = self.GTOs[a]
        B = self.GTOs[b]
        return overlap1s(A.alpha,A.type,A.center,B.alpha,B.type,B.center)

    def kineticInt(self,a,b):
        A = self.GTOs[a]
        B = self.GTOs[b]
        return kineticInt1s(A.alpha,A.type,A.center,B.alpha,B.type,B.center)
    
    def nucInt(self,a,b,R_C,Z):
        A = self.GTOs[a]
        B = self.GTOs[b]
        return nucInt1s(A.alpha,A.type,A.center,B.alpha,B.type,B.center,R_C,Z)

    def twoElecInt(self,a,b,c,d):
        A = self.GTOs[a]
        B = self.GTOs[b]
        C = self.GTOs[c]
        D = self.GTOs[d]
        return twoElecInt1s(A.alpha,A.type,A.center,B.alpha,B.type,B.center,
                                  C.alpha,C.type,C.center,D.alpha,D.type,D.center)


class RepGTO(Representation):
    """
    Representation class for the general GTO basis. 
    """
    def __init__(self,GTOs,Zs,nuclearPositions):
        self.GTOs = GTOs
        Representation.__init__(self,Zs,nuclearPositions,len(GTOs))

    def overLapInt(self,a,b):
        A = self.GTOs[a]
        B = self.GTOs[b]
        return overlap(A.alpha,A.type,A.center,B.alpha,B.type,B.center)
    
    def kineticInt(self,a,b):
        A = self.GTOs[a]
        B = self.GTOs[b]
        return kinetic(A.alpha,A.type,A.center,B.alpha,B.type,B.center)
    
    def nucInt(self,a,b,R_C,Z):
        A = self.GTOs[a]
        B = self.GTOs[b]
        return -Z*nuclear_attraction(A.alpha,A.type,A.center,B.alpha,B.type,B.center,R_C)
    
    def twoElecInt(self,a,b,c,d):
        A = self.GTOs[a]
        B = self.GTOs[b]
        C = self.GTOs[c]
        D = self.GTOs[d]
        #order of c and b's here is due to notational differences
        return electron_repulsion(A.alpha,A.type,A.center,C.alpha,C.type,C.center,
                                  B.alpha,B.type,B.center,D.alpha,D.type,D.center)


class RepCGTO(Representation):
    def __init__(self,GTOs,Zs,nuclearPositions):
        self.GTOs = GTOs
        Representation.__init__(self,Zs,nuclearPositions,len(GTOs))

        #should use RepGTO to get H matrix and then perform inner products to get contraction

    def overLapInt(self,a,b):  
        pass
    
    def kineticInt(self,a,b):
        pass
    
    def nucInt(self,a,b,R_C,Z):
        pass

    def twoElecInt(self,a,b,c,d):
        pass

