import sys
import os
import numpy as np
import matplotlib.pyplot as plt


sys.path.append(os.path.abspath("."))
from useful_code.numerov_method import *
from useful_code.verlet_method import *
"""
Should remove function .argmin() part and just work with index form
"""

def fKS(E,hatreePot):
    """
    f function for KS eq
    """
    return lambda r: -2*(E + 2/r - hatreePot(r))  

def uAsym(r):
    """
    Finds asymptotic solution
    """
    return np.sqrt(32)*r*np.exp(-2*r) #need to add correct normalisation

#is this the correct solution, does it need to be normalised?
#factor here doesn't seem to make a difference to intial conditions
#also coverges very quickly, doesn't seem to depend much on density
h = -0.01 #sign?
N = 999
acc = 0.001

r0 = 10
u0 = uAsym(r0)
r1 = r0 + h
u1 = uAsym(r1)

def getVs(uFunc):    #better name?
    Us = np.empty(N)
    Rs = np.empty(N) #use numpy arrays
    Us[0] = 0
    Rs[0] = 0
    Us[1] = 0.01
    Rs[1] = 0.01
    #should clean up here
    for i in range(2,N): #numerov method is not applicable here
        Us[i] = iterateVerlet(Us[i-1],Us[i-2],0.01,Rs[i-1],lambda U,r: -uFunc(r)**2/r if r !=0 else 0)
        Rs[i] = Rs[i-1] + 0.01

    #print(-1)
    #print(Us[-1])
    Us =  Us/abs(Us[-1]) #negative solutions here, problem?
    Vs = np.array(list(map(lambda r,U: U/r if r != 0 else 0,Rs,Us)))

    return Rs,Vs

def getVsFunc(rs,Vs):
    return lambda r: Vs[np.abs(rs-r).argmin()]

def getuFunc(rs,us):
    return lambda r: us[np.abs(rs-r).argmin()]


def findEigenstate(rs,us,EPrev):
    ERange = -np.arange(0.5,2,0.01)
    Rs,Vs = getVs(getuFunc(rs,us)) #there is a simpler way to structure this
    
    def findU0(E):
        rs0,us0 = runNumerov(r0,h,u0,u1,fKS(E,getVsFunc(Rs,Vs)),N-2) #better names
        if np.any(us0<0):
            return 10**10
        return us0[-1]
    epsilon = min(ERange,key=findU0) #rename

    rsNew,usNew = runNumerov(r0,h,u0,u1,fKS(epsilon,getVsFunc(Rs,Vs)),N-2) # should move to specifying with index?
    usNewNorm = usNew/np.sqrt(np.sum(abs(h)*usNew**2))
    RsNew,VsNew = getVs(getuFunc(rsNew,usNewNorm))
    ENew = 2*epsilon - np.sum(abs(h)*VsNew*usNewNorm**2)
    #print("-----")
    #print(epsilon)
    #print(np.sum(abs(h)*VsNew*usNewNorm**2))



    if abs(ENew - EPrev) < acc:
        return ENew
    else:
        return findEigenstate(rsNew,usNewNorm,ENew)


ERange = -np.arange(0.5,3,0.01)    

def findU0(E):
    rs0,us0 = runNumerov(r0,h,u0,u1,fKS(E,lambda r: 0),N-2)
    if np.any(us0[-1]<0):
        return 10**10
    return us0[-1] #geting some negative values

epsilon = min(ERange,key=findU0)
print(epsilon)
rs,us = runNumerov(r0,h,u0,u1,fKS(epsilon,lambda r: 0),N-2)
print(findEigenstate(rs,us/np.sqrt(np.sum(abs(h)*us**2)),epsilon))