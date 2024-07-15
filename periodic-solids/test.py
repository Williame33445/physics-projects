import numpy as np
from scipy import special as sp
import matplotlib.pyplot as plt

# constants
a = 6.822 #au
R = 0.5 #changes sign convention and increases the accuracy
vol = 3*(a**3)/4

l_max = 2*(2+1) +2
nLim = 2 

b1 = 2*np.pi*np.array([-1,1,1])/a
b2 = 2*np.pi*np.array([1,-1,1])/a
b3 = 2*np.pi*np.array([1,1,-1])/a
Ks =[l*b1 + m*b2 + n*b3 for l in range(-nLim,nLim+1) for m in range(-nLim,nLim+1) for n in range(-nLim,nLim+1)]

#generate A
A = np.empty([len(Ks),len(Ks)])
for i,K_i in enumerate(Ks):
    for j,K_j in enumerate(Ks):
        kDiff = np.linalg.norm(K_i-K_j)
        if kDiff == 0:
            A[i,j] = 1
        else:
            A[i,j] = -4*np.pi*(R**2)*sp.spherical_jn(1,kDiff*R)/(vol*kDiff)



#shall use hydrogen atom approximation as all other electrons have a shielding effect (this is rough)

def iterateRK4(t_n,y_n,h,f):
    """
    Iterates by RK4 method.
    """
    k_1 = f(t_n, y_n)
    k_2 = f(t_n + h/2, y_n + h*k_1/2)
    k_3 = f(t_n + h/2, y_n + h*k_2/2)
    k_4 = f(t_n + h, y_n + h*k_3)

    return y_n + h*(k_1 + 2*k_2 + 2*k_3 + k_4)/6

def intialCond(l,r_0):
    #from small r solution to radial schrodinger eq, remember its u not R
    return np.array([r_0**(l+1),(l+1)*r_0**l])

def V(r,Z):
    return -Z/r 

n = 10
def getfRadial(l,E,elecV,Z):
    def f(r,y):
        u = y[0]
        p = y[1] #p = du/dr
        f_u = p
        f_p = (l*(l+1)/(r**2) - 2*E + 2*(V(r,Z)+elecV))*u
        return np.array([f_u,f_p])
    return f
        
r_0 = 0.05
h = 0.05
lim = 0.1

def getRatio(l,E,Z,prev=0,elecV=0,tick=0):
    rs = [r_0]
    ys = [intialCond(l,r_0)]
    for i in range(int((R - r_0)/h)):
        rs.append(rs[-1]+h)
        ys.append(iterateRK4(rs[-1],ys[-1],h,getfRadial(l,E,elecV,Z)))

    new = ys[-1][0]
    if abs(prev - new) < lim or tick > 100:
#        if tick>100:
#            print("over")
        return ys[-1][1]/ys[-1][0] - 1/R
    
    elecV = sum([(h*y[0]**2)/rs[i] for i,y in enumerate(ys)])
    return getRatio(l,E,Z,new,elecV,tick+1)
def findDet(k,E,Z):
    #set up k dependent matricies 
    B = np.empty([len(Ks),len(Ks)])
    C = np.empty([len(Ks),len(Ks),len(Ks)])
    qs = [k + K for K in Ks]

    #find B and C
    for i in range(len(Ks)):
        for j in range(len(Ks)):
            qDot = np.dot(qs[i],qs[j])
            B[i,j] = A[i,j]*qDot/2
            qi = np.linalg.norm(qs[i])
            qj = np.linalg.norm(qs[j])
            for l in range(l_max): 
                if qi == 0 or qj == 0:
                    C[i,j,l] = 0
                else: 
                    C[i,j,l] = (2*l + 1)*2*np.pi*(R**2)*sp.eval_legendre(l,qDot/(qi*qj))*sp.spherical_jn(l,qi*R)*sp.spherical_jn(l,qj*R)/vol
    
    H = -E*A + B
    for l in range(l_max):
        H += C[:,:,l]*getRatio(l,E,Z)

    #this function is more approriate for large matrices
    detData = np.linalg.slogdet(H - E*np.identity(len(Ks))) #returns in (sign,absolute val
    det = detData[0]*detData[1]
    return det

EMesh  = np.arange(0,0.25,step=0.01) #np.arange(-0.02,0.25,step=0.01)
Z_d=13.201
print(findDet(2*np.pi*np.array([0.1,0,0])/a,0.1,Z_d))
#print(f"0.1:{EMesh[int(findZero([findDet(2*np.pi*np.array([0.1,0,0])/a,E,Z_d) for E in EMesh]))]}")
