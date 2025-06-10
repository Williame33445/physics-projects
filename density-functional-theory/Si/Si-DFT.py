#DFT program for local pseudopotential 

import numpy as np
from scipy.special import sph_harm

def sphericalHarmonic(l, m, k):
    phi = np.arctan(k[1]/k[0])
    theta = np.arccos(k[2]/np.linalg.norm(k))
    return sph_harm(m, l, phi, theta)

def GTHPseudopotential(zeta,Z_eff,V,C_1,C_2,r_s,r_p,h1s,h2s,h1p):

    def V_core(K): # 0 error for K=0
        KMag = np.linalg.norm(K)
        return -4*np.pi*Z_eff*np.exp(-0.5*(zeta*KMag)**2)/(V*KMag**2)
    
    def V_loc(K):
        KMag = np.linalg.norm(K)
        return (2*np.pi)*1.5*zeta**3*np.exp(-0.5*(zeta*KMag)**2)*(C_1 + C_2*(3 - (zeta*KMag)**2))/V
    
    def p1s(K):
        KMag = np.linalg.norm(K)
        return 4*r_s*(2*r_s/V)**0.5*np.pi**1.25*np.exp(-0.5*(KMag*r_s)**2)
    
    def p2s(K):
        KMag = np.linalg.norm(K)
        return 8*r_s(2*r_s/(15*V))**0.5*np.pi**1.25*np.exp(-0.5*(KMag*r_s)**2)*(3 - (KMag*r_s)**2)
    
    def p1p(K):
        KMag = np.linalg.norm(K)
        return 8*r_p**2*(r_p/3)**0.5*np.pi**1.25*np.exp(-0.5*(KMag*r_p)**2)*KMag/V
    
    def H_nonlocal(K1,K2):
        H = sphericalHarmonic(0,0,K1)*p1s(K1)*h1s*p1s(K2)*np.conjugate(sphericalHarmonic(0,0,K2))
        H += sphericalHarmonic(0,0,K1)*p2s(K1)*h2s*p2s(K2)*np.conjugate(sphericalHarmonic(0,0,K2))
        for i in range(-1,2):
            H -= sphericalHarmonic(1,i,K1)*p1p(K1)*h1p*p1p(K2)*np.conjugate(sphericalHarmonic(1,i,K2))
        return H
 
    def getH(K1,K2):
        return V_core(K1-K2) + V_loc(K1-K2) + H_nonlocal(K1,K2)
    
    return getH




    

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