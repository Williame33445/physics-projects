import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath("."))
from useful_code.numerov_method import *
from useful_code.verlet_method import *

def fKS(E,hatreePot):
    """
    f function for KS eq
    """
    return lambda r: -2*(E + 2/r - hatreePot(r))  

def uAsym(r):
    """
    Finds asymptotic solution
    """
    return np.sqrt(32)*r*np.exp(-2*r)

h = -0.01 
N = 999
acc = 0.001

r0 = 10
u0 = uAsym(r0)
r1 = r0 + h
u1 = uAsym(r1)

def getVHFunc(rs,us):
    uFunc = lambda r: us[np.abs(rs-r).argmin()]

    Us = np.empty(N)
    Rs = np.empty(N)
    Us[0] = 0
    Rs[0] = 0
    Us[1] = 0.01
    Rs[1] = 0.01

    for i in range(2,N):
        Us[i] = iterateVerlet(Us[i-1],Us[i-2],0.01,Rs[i-1],lambda U,r: -(uFunc(r)**2)/r)
        Rs[i] = Rs[i-1] + 0.01

    #there seems to be a negative problem here

    Us =  Us/Us[-1] #negative solutions here, problem?
    Vs = np.array(list(map(lambda r,U: U/r if r != 0 else 0,Rs,Us)))

    return lambda r: Vs[np.abs(Rs-r).argmin()]




def findEigenstate(rs,us,EPrev):
    print(EPrev)
    ERange = -np.arange(0.5,3,0.01)
    VsFunc = getVHFunc(rs,us) #may be able to hold in recursion
    
    def findU0(E):
        rs0,us0 = runNumerov(r0,h,u0,u1,fKS(E,VsFunc),N-2) #better names
        if np.any(us0<0): #removes unphysical solutions
            return 10**10
        return us0[-1]
    epsilon = min(ERange,key=findU0) 

    rsNew,usNew = runNumerov(r0,h,u0,u1,fKS(epsilon,VsFunc),N-2) # should move to specifying with index?
    usNewNorm = usNew/np.sqrt(np.sum(abs(h)*usNew**2))
    VsNewFunc = getVHFunc(rsNew,usNew)
    ENew = 2*epsilon - np.sum(abs(h)*VsNewFunc(rsNew)*(usNew)**2) 


    if abs(ENew - EPrev) < acc:
        return ENew
    else:
        return findEigenstate(rsNew,usNewNorm,ENew)


rs,us = runNumerov(r0,h,u0,u1,fKS(-1,lambda r: 1),N-2)

print(findEigenstate(rs,us/np.sqrt(np.sum(abs(h)*us**2)),-2.5))