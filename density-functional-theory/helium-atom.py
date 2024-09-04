import numpy as np
import matplotlib.pyplot as plt

#this program implements Thijssen's density function approach to the helium atom, it is currently broken

secMin = 0.0001 
def iterateSecant(x_0,x_1,f):
    """
    Applies secant method.
    """
    x_2 = (x_1*f(x_0) -x_0*f(x_1))/(f(x_0) - f(x_1))
    if abs(f(x_2)) < secMin:
        return x_2
    elif f(x_2) > 0:
        return iterateSecant(x_0,x_2,f)
    else:
        return iterateSecant(x_2,x_1,f)

def getf(Vs):
    """
    Finds us0 for sim.
    """
    def f(E):
        rs0,us0 = run(E,Vs)
        return us0[-1]
    return f



def uAsym(r):
    """
    Finds asymptotic solution
    """
    return r*np.exp(-r)

h = -0.01 
N = 500
acc = 0.001

r0 = 5 + h
u0 = uAsym(r0)
r1 = r0 + h
u1 = uAsym(r1)

def getVs(rs,us):
    """
    Finds Vs matrix
    """
    us = np.flip(us)
    Us = np.empty(N)
    Rs = np.empty(N)
    Us[0] = 0
    Rs[0] = 0
    Us[1] = 0.01
    Rs[1] = 0.01

    for i in range(2,N):
        Us[i] = 2*Us[i-1] - Us[i-2] - ((us[i-1]**2)/Rs[i-1])*(0.01**2)
        Rs[i] = Rs[i-1] + 0.01

    #negative solution problem occurs here

    Vs = np.array(list(map(lambda r,U: U/r if r != 0 else 0,Rs,Us))) #np way of doing this?

    return np.flip(Vs)

def run(E,Vs):
    """
    Solves for a given Vs
    """
    us = np.empty(N)
    rs = np.empty(N)
    us[0] = u0
    rs[0] = r0
    us[1] = u1
    rs[1] = r1

    for i in range(2,N):
        us[i] = 2*us[i-1] - us[i-2] - 2*(E + 2/rs[i-1] - Vs[i-1])*us[i-1]*(h**2)
        rs[i] = rs[i-1] + h
    usNorm = us/np.sqrt(np.sum(abs(h)*us**2))
    return rs, usNorm



def findEigenstate(rs,us,EPrev):
    """
    Recursive function that tries to deduce ground state.
    """
    Vs = getVs(rs,us)
    
    epsilon = iterateSecant(-2,-1,getf(Vs))

    rsNew,usNew = run(epsilon,Vs)
    VsNew = getVs(rsNew,usNew)
    ENew = 2*epsilon - np.sum(abs(h)*VsNew*(usNew**2))


    if abs(ENew - EPrev) < acc:
        return ENew
    return findEigenstate(rsNew,usNew,ENew)


rs,us = run(-2,np.zeros(N))
print(findEigenstate(rs,us,-2))