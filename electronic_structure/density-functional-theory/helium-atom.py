import sys
import os
import numpy as np

sys.path.append(os.path.abspath("."))
from useful_code.numerov_method import *

def fKS(E,hatreePot):
    """
    f function for KS eq
    """
    return lambda r: -2*(E + 2/r - hatreePot(r))

def fHartree(r,rs,us): #this is not the most efficient way of doing this
    """
    f function for Hartree eq
    """
    rMinIndex = np.abs(rs-r).argmin()
    return -(us[rMinIndex]**2)/rs[rMinIndex]

def uAsym(r):
    """
    Finds asymptotic solution
    """
    return np.sqrt(32)*r*np.exp(-2*r) #is this the correct solution, does it need to be normalised?

h = -0.1
N = 999
acc = 0.01

r0 = 100
u0 = uAsym(r0)
r1 = r0 + h
u1 = uAsym(r1)

def getVsFunc(rs,us):
    rs = np.flip(rs)
    us = np.flip(us)
    #print(us)
    rs,Us = runNumerov(0,abs(h),0,abs(h),lambda r: fHartree(r,rs,us),N)
    Us =  Us/Us[-1]
    def Vs(r):
        if np.all(r == 0):
            return 0
        rIndex = np.abs(rs-r).argmin()
        return Us[rIndex]/rs[rIndex]
    return Vs


def findEigenstate(rs,us,EPrev):
    print(EPrev)
    ERange = -np.arange(0.01,2,0.01)    

    def findU0(E):
        rs0,us0 = runNumerov(r0,h,u0,u1,fKS(E,getVsFunc(rs,us)),N)
        return us0[-1]
    epsilon = min(ERange,key=findU0) #rename

    rsNew,usNew = runNumerov(r0,h,u0,u1,fKS(epsilon,getVsFunc(rs,us)),N) # should move to specifying with index?
    usNewNorm = usNew/np.sqrt(np.sum(abs(h)*usNew**2))
    VsNew = getVsFunc(rsNew,usNewNorm)(rsNew)    
    ENew = 2*epsilon - np.sum(abs(h)*VsNew*usNewNorm**2)


    if abs(ENew - EPrev) < acc:
        #print(rsNew)
        #print(usNew)
         #need to simplify

        return ENew
    else:
        return findEigenstate(rsNew,usNewNorm,ENew)


def findUs(rs,us): #swapped order of arrays
    rs = np.flip(rs)
    us = np.flip(us)
    #print(us)
    rs,Us = runNumerov(0,abs(h),0,abs(h),lambda r: fHartree(r,rs,us),N)
    #print(Us[-1])
    return Us/Us[-1]   #alpha trick

rs,us = runNumerov(r0,h,u0,u1,fKS(-2,lambda r: 0),N)
print(findEigenstate(rs,us/np.sqrt(np.sum(abs(h)*us**2)),-2))