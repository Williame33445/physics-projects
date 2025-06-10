#DFT program 

import numpy as np
import GTHPseudopotential as GTH
from Atom import *

def u(r,q):
    return (2*r - q - 1)/(2*q)

class DFTSimulation:
    def __init__(self, ax, ay, az, KxMax, KyMax, KzMax, kxSampNum, kySampNum, kzSampNum, atoms, occuptationFunction):
        self.atoms = atoms
        self.ocurptationFunction = occuptationFunction

        ix = np.arange(-KxMax, KxMax+1)
        iy = np.arange(-KyMax, KyMax+1)
        iz = np.arange(-KzMax, KzMax+1)
        
        self.Ks = []
        for ix in range(-KxMax, KxMax+1):
            for iy in range(-KyMax, KyMax+1):
                for iz in range(-KzMax, KzMax+1):
                    if (ix/KxMax)**2 + (iy/KyMax)**2 + (iz/KzMax)**2 <= 1:
                        self.Ks.append(2*np.pi*np.array([ix/ax, iy/ay, iz/az]))

        self.getAtomicPotential()

        bx = 2*np.pi*np.array([1/ax,0,0])
        by = 2*np.pi*np.array([0,1/ay,0])
        bz = 2*np.pi*np.array([0,0,1/az])

        self.BZSample = []
        for rx in range(1, kxSampNum+1):
            for ry in range(1, kySampNum+1):
                for rz in range(1, kzSampNum+1):
                    self.BZSample.append(u(rx,kxSampNum)*bx + u(ry,kySampNum)*by + u(rz,kzSampNum)*bz)


    
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
            

        eigenvalues, eigenvectors = np.linalg(H)

        return eigenvalues, eigenvectors
    
    def getDensity(self, k):
        #the problem is how to calculate sum over Cs effeicently
        pass


HPseudoPotential = GTH.pseudopotential(0.2,1,125,-4.0663326,0.6778322,0,0,0,0,0)
atomList =  [Atom(np.array([0,0,0]),HPseudoPotential)]

DFTSimulation(5,5,5,1,1,1,atomList).solveForConstantDensity(np.array([0,0,0]))


#need to set up Si pseudopotential properly, also look at calculating density via BZ sampling