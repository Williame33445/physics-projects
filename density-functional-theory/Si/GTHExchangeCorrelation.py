import numpy as np

a1 = 0.4581652932831429
a2 = 2.217058676663745
a3 = 0.7405551735357053
a4 = 0.01968227878617998
b1 = 1.0
b2 = 4.504130959426697
b3 = 1.110667363742916
b4 = 0.02359291751427506

def Vxc(n):
    n = np.asarray(n)
    V = np.zeros_like(n, dtype=np.complex128)
    
    #takes nonzero elements of n
    mask = np.abs(n) > 1e-8
    n_nonzero = n[mask] 
    
    rs = (3 / (4 * np.pi * n_nonzero))**(1/3)

    den1 = -3*a1 -2*a2*rs -a3*rs**2
    num = b1*rs + b2*rs**2 + b3*rs**3 + b4*rs**4

    den21 = a1 + a2*rs + a3*rs**2 + a4 *rs**3
    den22 = b1*rs + 2*b2*rs**2 + 3*b3*rs**3 + 4*b4*rs**4
    
    #maps nonzero calculations to V arrayß
    V[mask] = (den1/num - den21*den22/num**2)/3

    return V