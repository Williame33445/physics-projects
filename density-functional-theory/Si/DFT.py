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

        KIxVals = np.arange(-KIxMax, KIxMax + 1)
        KIyVals = np.arange(-KIyMax, KIyMax + 1)
        KIzVals = np.arange(-KIzMax, KIzMax + 1)

        KIxGrid, KIyGrid, KIzGrid = np.meshgrid(KIxVals, KIyVals, KIzVals, indexing="ij")
        self.KIs = np.stack((KIxGrid, KIyGrid, KIzGrid), axis=-1).reshape(-1, 3)
        self.Ks = 2*np.pi*np.einsum("ij,kj",np.diag([1/ax, 1/ay, 1/az]), self.KIs).T

        self.getKsIndex = lambda K: (K[0] + KIxMax)*(2*KIxMax + 1)*(2*KIyMax + 1) + (K[1] + KIyMax)*(2*KIxMax + 1) + K[2] + KIzMax
        self.Ks0Index = self.getKsIndex(np.array([0, 0, 0]))


        self.getAtomicPotential()

        rxVals = np.arange(1, kxSamp + 1)
        ryVals = np.arange(1, kySamp + 1)
        rzVals = np.arange(1, kzSamp + 1)

        kxGrid, kyGrid, kzGrid = np.meshgrid(rxVals, ryVals, rzVals, indexing="ij")
        self.BZSamples = 2*np.pi*np.stack((u(kxGrid,kxSamp)/ax, u(kyGrid,kySamp)/ay, u(kzGrid,kzSamp)/az), axis=-1).reshape(-1, 3)
        
        

    
    def getAtomicPotential(self):
        self.V = np.zeros((len(self.Ks), len(self.Ks)), dtype=np.complex128)

        for i in range(len(self.Ks)):
            for j in range(i+1):
                for a in self.atoms:

                    self.V[i, j] += np.exp(-1j*np.dot(self.Ks[i] - self.Ks[j],a.position))*a.pseudopotential(self.Ks[i], self.Ks[j])
                self.V[j, i] = np.conjugate(self.V[i, j])

    
    def solveForConstantDensity(self,k): #can calculate kinetic part at start as well
        H = self.V.copy()

        for i in range(len(self.Ks)):
            H[i, i] = 0.5 * ((k[0] + self.Ks[i][0])**2 + (k[1] + self.Ks[i][1])**2 + (k[2] + self.Ks[i][2])**2)

        eigenvalues, eigenvectors = np.linalg.eig(H)

        return eigenvalues, eigenvectors
    
    def getDensity(self, KI):
        density = 0
        KIndexRel0 = self.getKsIndex(KI) - self.Ks0Index 
        for k_BZ in self.BZSamples:
            eigenvalues, eigenvectors = self.solveForConstantDensity(k_BZ)
            for i,e in enumerate(eigenvalues):
                occ = self.occuptationFunction(e)
                for KI_1 in self.KIs:
                    if np.all(np.abs(KI + KI_1) <= self.KIs[-1]): 
                        K1Index = self.getKsIndex(KI_1)
                        density += occ * np.conjugate(eigenvectors[i, KIndexRel0 + K1Index]) * eigenvectors[i, K1Index]/(len(self.BZSamples)*self.volume)
        return density




HPseudoPotential = GTH.pseudopotential(0.2,1,125,-4.0663326,0.6778322,0,0,0,0,0)
atomList =  [Atom(np.array([0,0,0]),HPseudoPotential)]

density = DFTSimulation(5,5,5,1,1,1,1,1,1,atomList,lambda e: 1 if e < 1 else 0).getDensity(np.array([0,0,0]))
print(density)