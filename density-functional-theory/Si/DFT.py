#DFT program 

import numpy as np
import GTHPseudopotential as GTHP
import GTHExchangeCorrelation as GTHEC
from Atom import *


def u(r,q):
    return (2*r - q - 1)/(2*q)

class DFTSimulation:
    def __init__(self, ax, ay, az, KIxMax, KIyMax, KIzMax, kxSamp, kySamp, kzSamp, atoms, occuptationFunction,Vxc):
        self.atoms = atoms
        self.occuptationFunction = occuptationFunction
        self.volume = ax*ay*az
        self.Vxc = Vxc

        #lattice sampling for plane wave basis
        self.KIxMax = KIxMax
        self.KIyMax = KIyMax
        self.KIzMax = KIzMax
        KIxs = np.arange(-KIxMax, KIxMax + 1)
        KIys = np.arange(-KIyMax, KIyMax + 1)
        KIzs = np.arange(-KIzMax, KIzMax + 1)
        self.KIxGrid, self.KIyGrid, self.KIzGrid = np.meshgrid(KIxs, KIys, KIzs, indexing="ij")

        #these lists give the K vectors associated with each index
        self.KIs = np.stack((self.KIxGrid, self.KIyGrid, self.KIzGrid), axis=-1).reshape(-1, 3)
        self.Ks = 2*np.pi*np.einsum("ij,kj",np.diag([1/ax, 1/ay, 1/az]), self.KIs).T

        # Function for finding method for find K index in eigenvector array
        self.getKsIndex = lambda KIx,KIy,KIz: (KIx + KIxMax)*len(KIxs)*len(KIys) + (KIy + KIyMax)*len(KIxs) + KIz + KIzMax
        self.Ks0Index = self.getKsIndex(0, 0, 0)

        #setting up BZ sampling
        rxs = np.arange(1, kxSamp + 1)
        rys = np.arange(1, kySamp + 1)
        rzs = np.arange(1, kzSamp + 1)
        kxGrid, kyGrid, kzGrid = np.meshgrid(rxs, rys, rzs, indexing="ij")
        self.BZSamples = 2*np.pi*np.stack((u(kxGrid,kxSamp)/ax, u(kyGrid,kySamp)/ay, u(kzGrid,kzSamp)/az), axis=-1).reshape(-1, 3)

        #intial density will be set to 0, this is just a simple calculation of density in this case
        self.getAtomicPotential()
        self.findFixedDensitySolution()

        self.nK = self.getKDensity(self.KIxGrid,self.KIyGrid,self.KIzGrid)

        VxcK = self.findKVxc(self.nK) #still need to associate with correct matrix elements, remeber k conventions
        #need to double the size of the sample region for the H
        #this just means setting up another K grid with double the dimentions in each direction and using this in the density calculation

        #need to find a way to make the different structures clearer        
        

    
    def getAtomicPotential(self):
        self.V = np.zeros((len(self.Ks), len(self.Ks)), dtype=np.complex128)
        for i in range(len(self.Ks)):
            for j in range(i+1):
                for a in self.atoms:
                    self.V[i, j] += np.exp(-1j*np.dot(self.Ks[i] - self.Ks[j],a.position))*a.pseudopotential(self.Ks[i], self.Ks[j])
                self.V[j, i] = np.conjugate(self.V[i, j])

    def getHamiltonian(self,k):
        H = self.V.copy()

        for i in range(len(self.Ks)):
            H[i, i] = 0.5 * ((k[0] + self.Ks[i][0])**2 + (k[1] + self.Ks[i][1])**2 + (k[2] + self.Ks[i][2])**2)

        return H
    
    def findFixedDensitySolution(self):
        self.epsilonss = []
        self.Css = []
        for k in self.BZSamples:
            H = self.getHamiltonian(k)
            eigenvalues, eigenvectors = np.linalg.eig(H)
            print(eigenvalues<1)
            self.epsilonss.append(eigenvalues)
            self.Css.append(eigenvectors)

    def getKDensity(self,KIxGrid, KIyGrid, KIzGrid): #self variables should be used here
        gridShape = self.KIxGrid.shape
        density = np.zeros(gridShape, dtype=np.complex128)

        for i in np.ndindex(gridShape):
            KIx = KIxGrid[i]
            KIy = KIyGrid[i]
            KIz = KIzGrid[i]
            KIndexRel0 = self.getKsIndex(KIx,KIy,KIz) - self.Ks0Index 
            for j,Cs in enumerate(self.Css):
                for k,e in enumerate(self.epsilonss[j]):
                    occ = self.occuptationFunction(e)
                    for KI_1 in self.KIs:
                        if abs(KIx + KI_1[0]) <= self.KIs[-1][0] and abs(KIy + KI_1[1]) <= self.KIs[-1][1] and abs(KIz + KI_1[2]) <= self.KIs[-1][2]:
                            K1Index = self.getKsIndex(*KI_1)
                            density[i] += occ * np.conjugate(Cs[k, KIndexRel0 + K1Index]) * Cs[k, K1Index]/(len(self.BZSamples)*self.volume)
        return density
    
    def findKVxc(self,nK): 
        ftKDensity = np.fft.fftn(nK)
        #may want to take a copy of grids to make the code clearer, try and find shorter names
        #phase introduced 
        ftPhaseFactor = np.exp(2*np.pi*1j*(self.KIxMax*(self.KIxGrid+self.KIxMax)/(1+2*self.KIxMax) + self.KIyMax*(self.KIyGrid+self.KIyMax)/(1+2*self.KIyMax) + self.KIzMax*(self.KIzGrid+self.KIzMax)/(1+2*self.KIzMax)))

        nP = ftPhaseFactor*ftKDensity
        #remember position conventions used here
        VxcP =  self.Vxc(nP) #remove values 
        VxcK = np.fft.ifftn(np.conjugate(ftPhaseFactor)*VxcP)
        return VxcK




HPseudoPotential = GTHP.pseudopotential(0.2,1,125,-4.0663326,0.6778322,0,0,0,0,0)
atomList =  [Atom(np.array([0,0,0]),HPseudoPotential)]

test = DFTSimulation(5,5,5,1,1,1,1,1,1,atomList,lambda e: 1 if e < 1 else 0,GTHEC.Vxc)
