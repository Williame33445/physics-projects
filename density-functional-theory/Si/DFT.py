#DFT program 

import numpy as np
import GTHPseudopotential as GTH
from Atom import *

def u(r,q):
    return (2*r - q - 1)/(2*q)

class DFTSimulation:
    def __init__(self, ax, ay, az, KxMax, KyMax, KzMax, kxSamp, kySamp, kzSamp, atoms, occuptationFunction):
        self.atoms = atoms
        self.occuptationFunction = occuptationFunction

        self.kxMax = KxMax
        self.kyMax = KyMax
        self.kzMax = KzMax


        KxVals = np.arange(-KxMax, KxMax + 1)
        KyVals = np.arange(-KyMax, KyMax + 1)
        KzVals = np.arange(-KzMax, KzMax + 1)

        self.getKIndex = lambda K: K[0]*(2*KxMax + 1)*(2*KyMax + 1) + K[1]*(2*KxMax + 1) + K[2] + (KzMax + KzMax + 1)


        KxGrid, KyGrid, KzGrid = np.meshgrid(KxVals, KyVals, KzVals, indexing="xy")
        self.Ks = np.stack((KxGrid, KyGrid, KzGrid), axis=-1).reshape(-1, 3)

        self.getAtomicPotential()

        self.bx = 2*np.pi*np.array([1/ax,0,0])
        self.by = 2*np.pi*np.array([0,1/ay,0])
        self.bz = 2*np.pi*np.array([0,0,1/az])

        kxVals = np.arange(1, kxSamp + 1)
        kyVals = np.arange(1, kySamp + 1)
        kzVals = np.arange(1, kzSamp + 1)

        kxGrid, kyGrid, kzGrid = np.meshgrid(kxVals, kyVals, kzVals, indexing="xy")
        self.BZSample = 2*np.pi*np.stack((u(kxGrid,kxSamp)/ax, u(kyGrid,kySamp)/ay, u(kzGrid,kzSamp)/az), axis=-1).reshape(-1, 3)
        

    
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
    
    def getDensity(self, K):
        density = 0
        for k_BZ in self.BZSample:
            eigenvalues, eigenvectors = self.solveForConstantDensity(k_BZ)
            for i,e in enumerate(eigenvalues):
                occ = self.occuptationFunction(e)
                for K_1 in self.Ks:
                    if np.all(np.abs(K + K_1) > self.KxMax + self.KyMax*self.by + self.kzMax*self.bz):
                        density += occ * np.conj(eigenvectors[i, self.getKIndex(K + K_1)]) * eigenvectors[i, self.getKIndex(K_1)]

        return density




HPseudoPotential = GTH.pseudopotential(0.2,1,125,-4.0663326,0.6778322,0,0,0,0,0)
atomList =  [Atom(np.array([0,0,0]),HPseudoPotential)]

eigenvalues, eigenvectors = DFTSimulation(5,5,5,1,1,1,1,1,1,atomList,lambda e: 1).solveForConstantDensity(np.array([0,0,0]))
print(eigenvalues)

#need to set up Si pseudopotential properly, also look at calculating density via BZ sampling