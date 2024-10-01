from abc import abstractmethod, ABC
import numpy as np 

from molecularIntegrals import GTOIntegrals 
from molecularIntegrals import GTO1sIntegrals


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
                integralij = self.kineticPlusNucInt(i,j)
                h[i,j] = integralij
                h[j,i] = integralij
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
    
    def normaliseElectronList(self,electrons):
        """
        Method that normalises a list electron objects.
        """
        for i,S in enumerate(electrons):
            electrons[i].ds = self.normalise(S.ds)
        return electrons
    
    def F(self,states1,states2):
        """
        Method that finds F^(+/-) matrix given C^(+/-) and C^(-/+) matricies.
        """
        sum1 = sum([np.einsum("prqs,r,s",self.twoElecInts, C, C) for C in states1])
        sum2 = sum([np.einsum("prsq,r,s",self.twoElecInts, C, C) for C in states1])
        sum3 = sum([np.einsum("prqs,r,s",self.twoElecInts, C, C) for C in states2])
        return self.h + sum1 - sum2 +sum3
    
    
    def findE(self,electrons):
        """
        Method that finds E given a list of electron objects.
        """
        return sum([C.energy + np.einsum("ij,i,j",self.h,C.ds,C.ds) for C in electrons])/2


    @abstractmethod
    def overLapInt(self,p,q):
        """
        Abstract method that finds the overlap integral elements.
        """
        pass

    @abstractmethod
    def kineticPlusNucInt(self,p,q):
        """
        Abstract method that finds the kinetic + nuclear integral elements.
        """
        pass

    @abstractmethod
    def twoElecInt(self,p,r,q,s):
        """
        Abstract method that finds the two electron integral elements.
        """
        pass


class RepGTO(Representation):
    """
    Representation class for primitive GTO basis. 
    """
    def __init__(self,GTOs,Zs,nuclearPositions):
        self.GTOs = GTOs
        Representation.__init__(self,Zs,nuclearPositions,len(GTOs))

    def overLapInt(self,a,b):
        A = self.GTOs[a]
        B = self.GTOs[b]
        return GTOIntegrals.overlap(A.alpha,A.type,A.center,B.alpha,B.type,B.center)
    
    def kineticPlusNucInt(self,a,b):
        A = self.GTOs[a]
        B = self.GTOs[b]
        kineticInt = GTOIntegrals.kinetic(A.alpha,A.type,A.center,B.alpha,B.type,B.center)
        nucInt = GTOIntegrals.nuc(A.alpha,A.type,A.center,B.alpha,B.type,B.center,self.nuclearPositions,self.Zs)
        return kineticInt + nucInt
    
    def twoElecInt(self,a,b,c,d):
        A = self.GTOs[a]
        B = self.GTOs[b]
        C = self.GTOs[c]
        D = self.GTOs[d]
        #order of c and b's here is due to notational differences
        return GTOIntegrals.electron_repulsion(A.alpha,A.type,A.center,C.alpha,C.type,C.center,
                                  B.alpha,B.type,B.center,D.alpha,D.type,D.center)


class RepCGTO(Representation):
    """
    Representation class for the contracted GTO basis.
    """
    def __init__(self,CGTOs,Zs,nuclearPositions):
        self.CGTOs = CGTOs

        primitiveGTOs = []
        self.CRepIndicies = []
        for CGTO in CGTOs:
            self.CRepIndicies.append(len(primitiveGTOs))
            primitiveGTOs += CGTO.primitives
        self.CRepIndicies.append(len(primitiveGTOs))

        self.primitiveRep = RepGTO(primitiveGTOs,Zs,nuclearPositions)

        Representation.__init__(self,Zs,nuclearPositions,len(CGTOs))
    
    def reduceMatrix(self,M,a,b):
        return M[self.CRepIndicies[a]:self.CRepIndicies[a+1],self.CRepIndicies[b]:self.CRepIndicies[b+1]]
    
    def reduce4Array(self,A,a,b,c,d):
        return A[self.CRepIndicies[a]:self.CRepIndicies[a+1],self.CRepIndicies[b]:self.CRepIndicies[b+1],
                 self.CRepIndicies[c]:self.CRepIndicies[c+1],self.CRepIndicies[d]:self.CRepIndicies[d+1]]

    def overLapInt(self,a,b):  
        return np.einsum("ij,i,j",self.reduceMatrix(self.primitiveRep.S,a,b),self.CGTOs[a].ds,self.CGTOs[b].ds) 
    
    def kineticPlusNucInt(self,a,b):
        return np.einsum("ij,i,j",self.reduceMatrix(self.primitiveRep.h,a,b),self.CGTOs[a].ds,self.CGTOs[b].ds) 
    
    def twoElecInt(self,a,b,c,d):
        return np.einsum("ijkl,i,j,k,l",self.reduce4Array(self.primitiveRep.twoElecInts,a,b,c,d),
                         self.CGTOs[a].ds,self.CGTOs[b].ds,self.CGTOs[c].ds,self.CGTOs[d].ds) 


class Rep1sGTO(Representation):
    """
    Representation class for only 1s GTOs.
    """
    def __init__(self,GTOs,Zs,nuclearPositions):
        self.GTOs = GTOs
        
        for GTO in GTOs:
            assert np.all(GTO.type == np.zeros(3))

        Representation.__init__(self,Zs,nuclearPositions,len(GTOs))

    def overLapInt(self,a,b):
        A = self.GTOs[a]
        B = self.GTOs[b]
        return GTO1sIntegrals.overlap1s(A.center,A.alpha,B.center,B.alpha)
    
    
    def kineticPlusNucInt(self,a,b):
        A = self.GTOs[a]
        B = self.GTOs[b]
        kineticInt = GTO1sIntegrals.kineticInt1s(A.center,A.alpha,B.center,B.alpha)
        nucInt = GTO1sIntegrals.nucInt1s(A.center,A.alpha,B.center,B.alpha,self.nuclearPositions,self.Zs)
        return kineticInt + nucInt
    
    def twoElecInt(self,a,b,c,d):
        A = self.GTOs[a]
        B = self.GTOs[b]
        C = self.GTOs[c]
        D = self.GTOs[d]
        return GTO1sIntegrals.twoElecInt2s(A.center,A.alpha,B.center,B.alpha,C.center,C.alpha,D.center,D.alpha)
