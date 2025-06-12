#DFT program 

import numpy as np
import GTHPseudopotential as GTH
from Atom import *

def u(r,q):
    return (2*r - q - 1)/(2*q)

class DFTSimulation:
    def __init__(self, ax, ay, az, KIxMax, KIyMax, KIzMax, kxSamp, kySamp, kzSamp, atoms, occuptationFunction):
        self.atoms = atoms
        self.occuptationFunction = occuptationFunction
        self.volume = ax*ay*az

        #lattice sampling for plane wave basis
        KIxVals = np.arange(-KIxMax, KIxMax + 1)
        KIyVals = np.arange(-KIyMax, KIyMax + 1)
        KIzVals = np.arange(-KIzMax, KIzMax + 1)
        KIxGrid, KIyGrid, KIzGrid = np.meshgrid(KIxVals, KIyVals, KIzVals, indexing="ij")
        self.KIs = np.stack((KIxGrid, KIyGrid, KIzGrid), axis=-1).reshape(-1, 3)
        self.Ks = 2*np.pi*np.einsum("ij,kj",np.diag([1/ax, 1/ay, 1/az]), self.KIs).T

        # Function for finding method for find K index in eigenvector array
        self.getKsIndex = lambda KI: (KI[0] + KIxMax)*len(KIxVals)*len(KIyVals) + (KI[1] + KIyMax)*len(KIxVals) + KI[2] + KIzMax
        self.Ks0Index = self.getKsIndex(np.array([0, 0, 0]))

        #setting up BZ sampling
        rxVals = np.arange(1, kxSamp + 1)
        ryVals = np.arange(1, kySamp + 1)
        rzVals = np.arange(1, kzSamp + 1)
        kxGrid, kyGrid, kzGrid = np.meshgrid(rxVals, ryVals, rzVals, indexing="ij")
        self.BZSamples = 2*np.pi*np.stack((u(kxGrid,kxSamp)/ax, u(kyGrid,kySamp)/ay, u(kzGrid,kzSamp)/az), axis=-1).reshape(-1, 3)

        #intial density will be set to 0, this is just a simple calculation of density in this case
        self.getAtomicPotential()
        self.findFixedDensitySolution()

        print(self.getDensity(np.array([0, 0, 0])))


        
        

    
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
            self.epsilonss.append(eigenvalues)
            self.Css.append(eigenvectors)

    def getDensity(self, KI):
        density = 0
        KIndexRel0 = self.getKsIndex(KI) - self.Ks0Index 
        for i,Cs in enumerate(self.Css):
            for i,e in enumerate(self.epsilonss[i]):
                occ = self.occuptationFunction(e)
                for KI_1 in self.KIs:
                    if np.all(np.abs(KI + KI_1) <= self.KIs[-1]): 
                        K1Index = self.getKsIndex(KI_1)
                        density += occ * np.conjugate(Cs[i, KIndexRel0 + K1Index]) * Cs[i, K1Index]/(len(self.BZSamples)*self.volume)
        return density




HPseudoPotential = GTH.pseudopotential(0.2,1,125,-4.0663326,0.6778322,0,0,0,0,0)
atomList =  [Atom(np.array([0,0,0]),HPseudoPotential)]

test = DFTSimulation(5,5,5,1,1,1,1,1,1,atomList,lambda e: 1 if e < 1 else 0)