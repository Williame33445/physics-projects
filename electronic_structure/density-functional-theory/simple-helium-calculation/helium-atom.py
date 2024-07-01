import numpy as np

"""
This program implements Thijssen's density function approach to the helium atom. This method is not 
particularly accurate. I have used integration methods for V_H instead of poisson equation methods.
"""

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

def findus0Func(Vs,V_x):
    """
    Finds us0 for sim.
    """
    def f(E):
        rs0,us0 = run(E,Vs,V_x) #change names
        return us0[-1]
    return f


h = -0.01 
N = 1000
acc = 0.001
lowerEpsilon = -2
upperEpsilon = -1

uAsym = lambda r: np.sqrt(32)*r*np.exp(-2*r)
r0 = 10 + h
u0 = uAsym(r0)
r1 = r0 + h
u1 = uAsym(r1)

def getV_H(rs,us):
    """
    Finds V_H, uses integration methods.
    """
    V_H = np.empty(N)
    V_H[-1] = np.sum((us**2)/rs)

    for i in range(N-1):
        V_H[i] = abs(h)*(np.sum((us[:i]**2)/rs[:i]) + np.sum(us[i:]**2)/rs[i])

    return 2*V_H  

def getV_x(rs,us):
    """
    Finds local density exchange potential.
    
    """
    return -(3*(us**2)/(2*(np.pi*rs)**2))**(1/3)


def run(epsilon,V_H,V_x):
    """
    Solves for a given V_H
    """
    us = np.empty(N)
    rs = np.empty(N)
    us[0] = u0
    rs[0] = r0
    us[1] = u1
    rs[1] = r1

    #this is the verlet method
    for i in range(2,N):
        us[i] = 2*us[i-1] - us[i-2] - 2*(epsilon + 2/rs[i-1] - V_H[i-1] - V_x[i-1])*us[i-1]*(h**2)
        rs[i] = rs[i-1] + h
    
    usNorm = us/np.sqrt(np.sum(abs(h)*us**2))
    return rs, usNorm



def findEigenstate(rs,us,EPrev):
    """
    Recursive function that tries to deduce ground state.
    """
    V_H = getV_H(rs,us)
    V_x = getV_x(rs,us)
    epsilon = iterateSecant(lowerEpsilon,upperEpsilon,findus0Func(V_H,V_x))

    rsNew,usNew = run(epsilon,V_H,V_x)
    V_HNew = getV_H(rsNew,usNew) 
    V_xNew = getV_x(rsNew,usNew)
    ENew = 2*epsilon - np.sum(abs(h)*V_HNew*(usNew**2)) + 0.5*np.sum(abs(h)*V_xNew*(usNew**2))

    if abs(ENew - EPrev) < acc:
        return ENew
    return findEigenstate(rsNew,usNew,ENew)


rs,us = run(-2,np.zeros(N),np.zeros(N))
print(findEigenstate(rs,us,-2))






