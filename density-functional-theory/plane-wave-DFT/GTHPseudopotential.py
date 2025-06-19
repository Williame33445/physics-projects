
import numpy as np
from scipy.special import sph_harm
    

def sphericalHarmonic(l, m, k):
    kMag = np.linalg.norm(k)
    if kMag == 0:
        return sph_harm(m, l, 0, 0)
    
    theta = np.real_if_close(np.arccos(k[2]/kMag))

    if k[0] == 0:
        phi = np.pi/2
    else:
        phi = np.real_if_close(np.arctan(k[1]/k[0]))

    return sph_harm(m, l, phi, theta)

def pseudopotential(zeta,Z_eff,V,C_1,C_2,r_s,r_p,h1s,h2s,h1p):

    def V_core(K): # 0 error for K=0
        KMag = np.linalg.norm(K)
        if KMag == 0:
            return 0
        return -4*np.pi*Z_eff*np.exp(-0.5*(zeta*KMag)**2)/(V*KMag**2)
    
    def V_loc(K):
        KMag = np.linalg.norm(K)
        return (2*np.pi)**1.5*zeta**3*np.exp(-0.5*(zeta*KMag)**2)*(C_1 + C_2*(3 - (zeta*KMag)**2))/V
    
    def p1s(K):
        KMag = np.linalg.norm(K)
        return 4*r_s*(2*r_s/V)**0.5*np.pi**1.25*np.exp(-0.5*(KMag*r_s)**2)
    
    def p2s(K):
        KMag = np.linalg.norm(K)
        return 8*r_s*(2*r_s/(15*V))**0.5*np.pi**1.25*np.exp(-0.5*(KMag*r_s)**2)*(3 - (KMag*r_s)**2)
    
    def p1p(K):
        KMag = np.linalg.norm(K)
        return 8*r_p**2*(r_p/(3*V))**0.5*np.pi**1.25*np.exp(-0.5*(KMag*r_p)**2)*KMag
    
    def H_nonlocal(K1,K2):
        H = sphericalHarmonic(0,0,K1)*p1s(K1)*h1s*p1s(K2)*np.conjugate(sphericalHarmonic(0,0,K2))
        H += sphericalHarmonic(0,0,K1)*p2s(K1)*h2s*p2s(K2)*np.conjugate(sphericalHarmonic(0,0,K2))
        for i in range(-1,2):
            H -= sphericalHarmonic(1,i,K1)*p1p(K1)*h1p*p1p(K2)*np.conjugate(sphericalHarmonic(1,i,K2))
        return H
 
    def getH(K1,K2):
        return V_core(K1-K2) + V_loc(K1-K2) + H_nonlocal(K1,K2)
    
    return getH
