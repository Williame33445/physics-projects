import sys
import os
import numpy as np

sys.path.append(os.path.abspath("."))
from useful_code.numerov_method import *

def f(E,hatreePot):
    return lambda r: -2*(E + 2/r - hatreePot(r))

def g(r,rs,us): #need better names
    rIndex = np.abs(rs-r).argmin()
    return -(us[rIndex]**2)/rs[rIndex]

def uAsym(r):
    return r*np.exp(-r) 

h = -0.1
N = 999
acc = 0.1

r0 = 100
u0 = uAsym(r0)
r1 = r0 + h
u1 = uAsym(r1)

def getVsFunc(rs,us):
    Us = findUs(rs,us) #not the most efficient way of doing this
    def Vs(r):
        if np.all(r == 0):
            return 0
        rIndex = np.abs(Us-r).argmin()
        return Us[rIndex]/np.flip(rs)[rIndex]
    return Vs


def findEigenstate(rs,us,EPrev):
    ERange = -np.arange(0.01,2,0.01)    

    def findU0(E):
        rs0,us0 = runNumerov(r0,h,u0,u1,f(E,getVsFunc(rs,us)),N)
        return us0[-1]
    E = min(ERange,key=findU0) #rename

    rsNew,usNew = runNumerov(r0,h,u0,u1,f(E,getVsFunc(rs,us)),N) # should move to specifying with index?
    usNewNorm = usNew/np.sqrt(np.sum(abs(h)*usNew**2))
    #print(usNewNorm)

    if abs(E - EPrev) < acc:
        #print(rsNew)
        #print(usNew)
        VsNew = np.vectorize(getVsFunc(rsNew,usNewNorm))(np.flip(rsNew))
        print(VsNew)
        print(E)
        
        Eco = 2*E - np.sum(abs(h)*VsNew*usNewNorm**2) #need to simplify
        print(2*E-Eco)

        return Eco
    else:
        return findEigenstate(rsNew,usNewNorm,E)


def findUs(rs,us): #swapped order of arrays
    rs = np.flip(rs)
    us = np.flip(us)
    #print(us)
    rs,Us = runNumerov(0,abs(h),0,1,lambda r: g(r,rs,us),N) #alpha value?
    return Us    

rs,us = runNumerov(r0,h,u0,u1,f(-2,lambda r: 0),N)
print(findEigenstate(rs,us/np.sqrt(np.sum(abs(h)*us**2)),-2))