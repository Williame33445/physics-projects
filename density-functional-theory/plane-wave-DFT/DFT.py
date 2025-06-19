#DFT program, only considers systems with an even number of electrons in the supercell.

import numpy as np
from Atom import *
from scipy.fft import fftn, ifftn


def u(r,q):
    return (2*r - q - 1)/(2*q)

class DFTSimulation:
    def __init__(self, ax, ay, az, KIxMax, KIyMax, KIzMax, kxSamp, kySamp, kzSamp, atoms, N,Vxc, tol):
        self.atoms = atoms
        self.N = N
        self.volume = ax*ay*az
        self.Vxc = Vxc
        self.ax = ax
        self.ay = ay
        self.az = az

        #lattice sampling for plane wave basis
        self.KHIxMax = KIxMax
        self.KHIyMax = KIyMax
        self.KHIzMax = KIzMax
        KHIxs = np.arange(-KIxMax, KIxMax + 1)
        KHIys = np.arange(-KIyMax, KIyMax + 1)
        KHIzs = np.arange(-KIzMax, KIzMax + 1)
        self.KHIxGrid, self.KHIyGrid, self.KHIzGrid = np.meshgrid(KHIxs, KHIys, KHIzs, indexing="ij")
        self.KnIxMax = 2*KIxMax
        self.KnIyMax = 2*KIyMax
        self.KnIzMax = 2*KIzMax
        KnIxs = np.arange(-2*KIxMax, 2*KIxMax + 1)
        KnIys = np.arange(-2*KIyMax, 2*KIyMax + 1)
        KnIzs = np.arange(-2*KIzMax, 2*KIzMax + 1)
        self.KnIxGrid, self.KnIyGrid, self.KnIzGrid = np.meshgrid(KnIxs, KnIys, KnIzs, indexing="ij")

        #these lists give the K vectors associated with each index
        self.KHIs = np.stack((self.KHIxGrid, self.KHIyGrid, self.KHIzGrid), axis=-1).reshape(-1, 3)
        self.KHs = 2*np.pi*np.einsum("ij,kj",np.diag([1/ax, 1/ay, 1/az]), self.KHIs).T

        # Function for finding method for find K index in eigenvector array
        self.getKHsIndex = lambda KIx,KIy,KIz: (KIx + KIxMax)*len(KHIxs)*len(KHIys) + (KIy + KIyMax)*len(KHIxs) + KIz + KIzMax
        self.KHs0Index = self.getKHsIndex(0, 0, 0)

        #setting up BZ sampling
        rxs = np.arange(1, kxSamp + 1)
        rys = np.arange(1, kySamp + 1)
        rzs = np.arange(1, kzSamp + 1)
        kxGrid, kyGrid, kzGrid = np.meshgrid(rxs, rys, rzs, indexing="ij")
        self.BZSamples = 2*np.pi*np.stack((u(kxGrid,kxSamp)/ax, u(kyGrid,kySamp)/ay, u(kzGrid,kzSamp)/az), axis=-1).reshape(-1, 3)


        self.HKins = [] 
        for k in self.BZSamples:
            self.HKins.append(self.getHKin(k))


        self.getAtomicPotential()
        self.nK = np.full((len(KnIxs),len(KnIys),len(KnIzs)), 1e-4, dtype=np.float64)

        #at each stage of the iterative cycle
        for iter in range(1000):
            VxcVha = self.getKVxcVha()
            self.Veff = VxcVha + self.Va 
            self.findFixedDensitySolution()

            newDensity = self.getNewDensity()
            if np.max(np.abs(newDensity - self.nK)) < tol and iter > 0:
                print("Converged")
                print(np.array(self.epsilonss))
                break
            else:
                print(np.max(np.abs(newDensity - self.nK)))
                mixingParameter = 0.9 if iter < 10 else 0.6
                self.nK = mixingParameter*newDensity + (1-mixingParameter)*self.nK

        

    
    def getAtomicPotential(self):
        self.Va = np.zeros((len(self.KHs), len(self.KHs)), dtype=np.complex128)
        for i in range(len(self.KHs)):
            for a in self.atoms:
                self.Va[i, i] += a.pseudopotential(self.KHs[i], self.KHs[i])
            for j in range(i):
                for a in self.atoms:
                    self.Va[i, j] += np.exp(-1j*np.dot(self.KHs[i] - self.KHs[j],a.position))*a.pseudopotential(self.KHs[i], self.KHs[j])
                self.Va[j, i] = np.conjugate(self.Va[i, j])

    def getHKin(self,k):
        HKin = np.zeros((len(self.KHs),len(self.KHs)),dtype=np.complex128)

        for i in range(len(self.KHs)): 
            HKin[i, i] = 0.5 * ((k[0] + self.KHs[i][0])**2 + (k[1] + self.KHs[i][1])**2 + (k[2] + self.KHs[i][2])**2)
        
        return HKin
    
    def findFixedDensitySolution(self):
        self.epsilonss = []
        self.Css = []
        for i in range(len(self.BZSamples)):
            H = self.HKins[i] + self.Veff
            eigenvalues, eigenvectors = np.linalg.eigh(H)
            self.epsilonss.append(eigenvalues)
            self.Css.append(eigenvectors)

    def getNewDensity(self): 
        newDensity = np.zeros_like(self.nK, dtype=np.complex128)
        epsilons = np.array(self.epsilonss).flatten()
        fermiEnergy = epsilons[epsilons.argsort()][int(self.N/2)-1]
        for i in np.ndindex(newDensity.shape):
            KIx = self.KnIxGrid[i]
            KIy = self.KnIyGrid[i]
            KIz = self.KnIzGrid[i]
            KIndexRel0 = self.getKHsIndex(KIx,KIy,KIz) - self.KHs0Index 
            for j,Cs in enumerate(self.Css):
                for k,e in enumerate(self.epsilonss[j]): 
                    for KI_1 in self.KHIs:
                        if -self.KHIxMax <= KIx + KI_1[0] <= self.KHIxMax and -self.KHIyMax <= KIy + KI_1[1] <= self.KHIyMax and -self.KHIzMax <= KIz + KI_1[2] <= self.KHIzMax:
                            K1Index = self.getKHsIndex(*KI_1)
                            newDensity[i] += 2*np.conjugate(Cs[KIndexRel0 + K1Index, k]) * Cs[K1Index, k]/(len(self.BZSamples)*self.volume) if e < fermiEnergy else 0
        return newDensity
    
    def getKVxcVha(self):
        pf1 = np.exp(2*np.pi*1j*(self.KnIxMax*(self.KnIxGrid + self.KnIxMax)/(1+2*self.KnIxMax) + self.KnIyMax*(self.KnIyGrid + self.KnIyMax)/(1+2*self.KnIyMax) + self.KnIzMax*(self.KnIzGrid + self.KnIzMax)/(1+2*self.KnIzMax)))
        pf2 = np.exp(2*np.pi*1j*(self.KnIxMax*self.KnIxGrid/(1+2*self.KnIxMax) + self.KnIyMax*self.KnIyGrid/(1+2*self.KnIyMax) + self.KnIzMax*self.KnIzGrid/(1+2*self.KnIzMax)))

        nP = pf2*fftn(pf1*self.nK)

        VxcP =  self.Vxc(nP) 
        VxcK = np.conjugate(pf1)*ifftn(np.conjugate(pf2)*VxcP) 
        
        VxcVha = np.zeros((len(self.KHs), len(self.KHs)), dtype=np.complex128)


        for i in range(len(self.KHs)):
            KIi = self.KHIs[i]
            VxcVha[i,i] = VxcK[2*self.KHIxMax, 2*self.KHIyMax, 2*self.KHIzMax]
            for j in range(i):
                KIj = self.KHIs[j]
                Vxc = VxcK[KIi[0]-KIj[0]+ 2*self.KHIxMax, KIi[1]-KIj[1] + 2*self.KHIyMax, KIi[2] - KIj[2] + 2*self.KHIzMax]
                Vha =  4*np.pi*self.nK[KIi[0]-KIj[0]+ 2*self.KHIxMax, KIi[1]-KIj[1]+ 2*self.KHIyMax, KIi[2] - KIj[2]+ 2*self.KHIzMax]/np.dot(self.KHs[i]-self.KHs[j],self.KHs[i]-self.KHs[j])
                VxcVha[i,j] = Vxc + Vha
                VxcVha[j,i] = np.conjugate(VxcVha[i,j])

        return VxcVha
    
    def findnP(self,xG,yG,zG):
        nP = np.zeros(xG.shape)
        for i in range(xG.shape[0]):
            for j in range(xG.shape[1]):
                for k in range(xG.shape[2]):
                    nP[i,j,k] = np.sum(self.nK*np.exp(-2*np.pi*1j*(self.KnIxGrid*xG[i,j,k]/self.ax + self.KnIyGrid*yG[i,j,k]/self.ay + self.KnIzGrid*zG[i,j,k]/self.az)))
        return nP
