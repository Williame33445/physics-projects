#DFT program 

import numpy as np
import GTHPseudopotential as GTHP
import GTHExchangeCorrelation as GTHEC
from Atom import *
import matplotlib.pyplot as plt



def u(r,q):
    return (2*r - q - 1)/(2*q)

class DFTSimulation:
    def __init__(self, ax, ay, az, KIxMax, KIyMax, KIzMax, kxSamp, kySamp, kzSamp, atoms, occuptationFunction,Vxc):
        self.atoms = atoms
        self.occuptationFunction = occuptationFunction
        self.volume = ax*ay*az
        self.Vxc = Vxc

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
        print(self.BZSamples)

        #intial density will be set to 0, this is just a simple calculation of density in this case
        self.getAtomicPotential()
        self.nK = np.zeros((len(KnIxs),len(KnIys),len(KnIzs)), dtype=np.float64)

        #at each stage of the iterative cycle
        for iter in range(100):
            VxcVha = self.getKVxcVha()
            self.Veff = VxcVha + self.V #should get kinetic part saved in here too
            self.findFixedDensitySolution()

            newDensity = self.getNewDensity()
            if np.max(np.abs(newDensity - self.nK)) < 1e-3:
                print("Converged")
                print(np.array(self.epsilonss))
                break
            else:
                print(np.max(np.abs(newDensity - self.nK)))
                self.nK = newDensity

        

    
    def getAtomicPotential(self):
        self.V = np.zeros((len(self.KHs), len(self.KHs)), dtype=np.complex128)
        for i in range(len(self.KHs)):
            for j in range(i+1):
                for a in self.atoms: #local and nonlocal part needs to be applied differently
                    self.V[i, j] += np.exp(-1j*np.dot(self.KHs[i] - self.KHs[j],a.position))*a.pseudopotential(self.KHs[i], self.KHs[j])
                self.V[j, i] = np.conjugate(self.V[i, j])

    def getHamiltonian(self,k):
        H = self.Veff.copy()

        for i in range(len(self.KHs)): #could make this more efficient by saving the kinetic part initally
            H[i, i] = 0.5 * ((k[0] + self.KHs[i][0])**2 + (k[1] + self.KHs[i][1])**2 + (k[2] + self.KHs[i][2])**2)
        
        return H
    
    def findFixedDensitySolution(self):
        self.epsilonss = []
        self.Css = []
        for k in self.BZSamples:
            H = self.getHamiltonian(k)
            eigenvalues, eigenvectors = np.linalg.eigh(H)
            #sorts eigenvalues from smallest to largest
            idx = eigenvalues.argsort()
            self.epsilonss.append(eigenvalues[idx])
            self.Css.append(eigenvectors[idx])

    def getNewDensity(self): 
        newDensity = np.zeros_like(self.nK, dtype=np.complex128)
        for i in np.ndindex(newDensity.shape):
            KIx = self.KnIxGrid[i]
            KIy = self.KnIyGrid[i]
            KIz = self.KnIzGrid[i]
            KIndexRel0 = self.getKHsIndex(KIx,KIy,KIz) - self.KHs0Index 
            Css = [self.Css[0]]
            for j,Cs in enumerate(Css):
                for k,e in enumerate(self.epsilonss[j]):
                    occ = self.occuptationFunction(e)
                    for KI_1 in self.KHIs:
                        if 0 <= KIx + KI_1[0] < 2*self.KHIxMax and 0 <= KIy + KI_1[1] < 2*self.KHIyMax and 0 <= KIz + KI_1[2] < 2*self.KHIzMax:
                            K1Index = self.getKHsIndex(*KI_1)
                            newDensity[i] += occ * np.conjugate(Cs[k, KIndexRel0 + K1Index]) * Cs[k, K1Index]/(len(self.BZSamples)*self.volume)
        print(np.all(0==newDensity-newDensity.T.conj())) #not hermitian from the first iteration (need to look at this)
        return newDensity
    
    def getKVxcVha(self):#the error is in here
        pf1 = np.exp(2*np.pi*1j*(self.KnIxMax*(self.KnIxGrid + self.KnIxMax)/(1+2*self.KnIxMax) + self.KnIyMax*(self.KnIyGrid + self.KnIyMax)/(1+2*self.KnIyMax) + self.KnIzMax*(self.KnIzGrid + self.KnIzMax)/(1+2*self.KnIzMax)))
        pf2 = np.exp(2*np.pi*1j*(self.KnIxMax*self.KnIxGrid/(1+2*self.KnIxMax) + self.KnIyMax*self.KnIyGrid/(1+2*self.KnIyMax) + self.KnIzMax*self.KnIzGrid/(1+2*self.KnIzMax)))

        #print(np.all(0==self.nK - self.nK.T))
        nP = pf2*np.fft.fftn(pf1*self.nK)
        self.nP = nP
        #remember position conventions used here
        VxcP =  self.Vxc(nP) #remove values 
        VxcK = np.conjugate(pf1)*np.fft.ifftn(np.conjugate(pf2)*VxcP) #clearer name
        
        VxcVha = np.zeros((len(self.KHs), len(self.KHs)), dtype=np.complex128)

        for i in range(len(self.KHs)):
            KIi = self.KHIs[i]
            for j in range(i+1):
                KIj = self.KHIs[j]
                Vxc = VxcK[KIi[0]-KIj[0]+ self.KHIxMax, KIi[1]-KIj[1] + self.KHIyMax, KIi[2] - KIj[2] + self.KHIzMax]
                if  not np.all(self.KHs[i] == self.KHs[j]):
                    Vha =  4*np.pi*self.nK[KIi[0]-KIj[0]+ self.KHIxMax, KIi[1]-KIj[1]+ self.KHIyMax, KIi[2] - KIj[2]+ self.KHIzMax]/np.dot(self.KHs[i]-self.KHs[j],self.KHs[i]-self.KHs[j])
                else:
                    Vha = 0


                VxcVha[i,j] = Vxc + Vha
                VxcVha[j,i] = np.conjugate(VxcVha[i,j])

        return VxcVha



ax = 5
ay = 5
az = 5
HPseudoPotential = GTHP.pseudopotential(0.2,1,ax*ay*az,-4.0663326,0.6778322,0,0,0,0,0)
atomList =  [Atom(np.array([2.5,2.5,2.5]),HPseudoPotential)]

test = DFTSimulation(ax,ay,az,1,1,1,1,1,1,atomList,lambda e: 1 if e < 0 else 0,GTHEC.Vxc)
#need to work out what is going wrong here, why is the 0 density part wrong
#this is not even to do with the convergence, its a problem with the intial part

print(test.nP) # complex density is a real issue

#need to redo occupations, rember that bands can overlap, my approach here assumes that this doesn't occur