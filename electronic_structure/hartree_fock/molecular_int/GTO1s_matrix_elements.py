import numpy as np
from scipy import special as sp
"""
This py file contains functions required to calculate matrix elements for symmetric GTO basis.
"""

def F_0(t):
    if t == 0:
        return 1 
    else:
        return np.sqrt(np.pi)*sp.erf(t)/(2*t)

def overlap1s(R_a,alpha_a,R_b,alpha_b):
    term1 = (np.pi/(alpha_a + alpha_b))**1.5
    term2 = np.exp(-(alpha_a*alpha_b)*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b))
    return term1*term2

def kineticInt1s(R_a,alpha_a,R_b,alpha_b):
    term1 = alpha_a*alpha_b/(alpha_a + alpha_b)
    term2 = 6 - 4*alpha_a*alpha_b*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b)
    term3 = (np.pi/(alpha_a + alpha_b))**1.5
    term4 = np.exp(-(alpha_a*alpha_b)*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b))

    return 0.5*term1*term2*term3*term4

def nucInt1s(R_a,alpha_a,R_b,alpha_b,R_C,Z):
    R_P = (alpha_a*R_a + alpha_b*R_b)/(alpha_a + alpha_b)
    term1 = -2*np.pi*Z/(alpha_a + alpha_b)
    term2 = np.exp(-(alpha_a*alpha_b)*(np.linalg.norm(R_a - R_b)**2)/(alpha_a + alpha_b))

    t = ((alpha_a + alpha_b)*(np.linalg.norm(R_P - R_C)**2))**0.5
    term3 = F_0(t)

    return term1*term2*term3

def twoElecInt2s(R_a,alpha_a,R_b,alpha_b,R_c,alpha_c,R_d,alpha_d):
    R_P = (alpha_a*R_a + alpha_c*R_c)/(alpha_a + alpha_c)
    R_Q = (alpha_b*R_b + alpha_d*R_d)/(alpha_b + alpha_d)

    term1 = 2*(np.pi**2.5)/((alpha_a + alpha_c)*(alpha_b + alpha_d)*((alpha_a + alpha_b + alpha_c + alpha_d)**0.5))
    term2 = np.exp(-(alpha_a*alpha_c*(np.linalg.norm(R_a - R_c)**2)/(alpha_a + alpha_c)) - (alpha_b*alpha_d*(np.linalg.norm(R_b - R_d)**2)/(alpha_b + alpha_d)))
    
    t = np.linalg.norm(R_P - R_Q)*((alpha_a + alpha_c)*(alpha_b + alpha_d)/(alpha_a + alpha_b + alpha_c + alpha_d))**0.5
    term3 = F_0(t)

    return term1*term2*term3
