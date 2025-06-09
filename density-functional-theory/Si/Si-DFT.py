#DFT program for local pseudopotential 

import numpy as np

def SiPseudopotential(k):
    #from Thijssen
    if np.isclose(k, np.sqrt(3)):
        return -0.1121
    elif np.isclose(k, np.sqrt(8)):
        return 0.0276
    elif np.isclose(k, np.sqrt(11)):
        return 0.0362
    else:
        return 0.0
    

class DFTSimulation:
    def __init__(self, Lx, Ly, Lz, KxMax, KyMax, KzMax,pseudopotential):
        self.pseudopotential = pseudopotential

        ix = np.arange(-KxMax, KxMax+1)
        iy = np.arange(-KyMax, KyMax+1)
        iz = np.arange(-KzMax, KzMax+1)
        
        self.Ks = []
        for ix in range(-KxMax, KxMax+1):
            for iy in range(-KyMax, KyMax+1):
                for iz in range(-KzMax, KzMax+1):
                    if (ix/KxMax)**2 + (iy/KyMax)**2 + (iz/KzMax)**2 <= 1:
                        self.Ks.append(2*np.pi*np.array([ix/Lx, iy/Ly, iz/Lz]))
    
    def solveForConstantDensity(self,k):
        M = np.zeros((len(self.Ks),len(self.Ks)), dtype=np.complex128)

        for i,K in enumerate(self.Ks):
            M[i, i] = 0.5 * ((k[0] + K[0])**2 + (k[1] + K[1])**2 + (k[2] + K[2])**2)

        print(len(np.linalg.eigvals(M)))


    

DFTSimulation(1,1,1,1,1,1,1).solveForConstantDensity(np.array([0,0,0]))

#need to set up Si pseudopotential properly, also look at calculating density via BZ sampling