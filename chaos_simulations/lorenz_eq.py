import sys
import os

sys.path.append(os.path.abspath("."))
from useful_code.RK4_methods import *

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#constants
sigma = 10
b = 8/3
r = 28

#defition of Lorenz eq
def f(t,R):
    x = R[0]
    y = R[1]
    z = R[2]
    f_x = sigma*(y - x)
    f_y = r*x - y - x*z
    f_z = x*y - b*z
    return np.array([f_x,f_y,f_z])

#simulation parameter
h = 0.001
N = 10**6

#simulate function
def simulate(RStates):
    t = 0
    for x in range(N):
        RStates.append(iterateRK4(t,RStates[-1],h,f))
    return RStates

def transformToXYZ(lst):
    return [[s[c] for s in lst] for c in range(3)]


#Deviation of intially close trajectories, try and get liapunov coefficent (need more data)

#delta = 10**(-10) 
#state1 = simulate([np.array([0.1,0,0])])
#state2 = simulate([np.array([0.1+delta,0,0])])
#deviationList = [np.linalg.norm(state1[i]-state2[i]) for i in range(len(state1))]
#timeArray =  np.linspace(0,1,num=N+1) ,need to scale to get actual predictions

#plt.plot(timeArray,deviationList)
#plt.show()




#predictive chaos structure (done, need to add some axes and clean up)

#niave way to find local maxima
def maxima(i,lst):
    if i == 0 or i == 1 or i == len(lst) - 1 or i == len(lst):
        return False
    return (lst[i-1]<lst[i]) and (lst[i]>lst[i+1])

#s1 = transformToXYZ(simulate([np.array([0.1,0,0])]))[2]
#maximaList = [m for i,m in enumerate(s1) if maxima(i,s1)]

#z_n = maximaList[:-1]
#z_n1 = maximaList[1:]
#plt.scatter(z_n,z_n1)
#plt.show()



# example simulation and graph
#states = simulate([np.array([0.1,0,0])])

#ax = plt.axes(projection='3d')
#ax.plot3D(states[0], states[1], states[2], 'blue',linewidth = '.5')
#plt.show()

# Use this to generate a large number of simulations (via GPUs?) then try and get dim via averages
num = 2
evenList = []
for x in np.linspace(-.5,.5,num):
    for y in np.linspace(-.5,.5,num):
        for z in np.linspace(-.5,.5,num):
            evenList.append(np.array([x,y,z]))

statesLst = []
for l in evenList:
    statesLst += simulate([l])

print(statesLst)
def findN(R,lst):
    n = 0 #should write this more elgantly
    for pos in lst:
        if np.linalg.norm(pos)<R:
            n += 1
    return n

R,N = np.linspace(0,20,40),np.linspace(0,20,40) #naming here is bad
for i,r in enumerate(R):
    N[i] = findN(r,statesLst)

plt.plot(R ,np.log(N))
plt.show()

#use method of least squares to deduce dimention
m,c = np.linalg.lstsq(np.vstack([R,np.ones(len(R))]).T,np.log(N),rcond=None)[0]
print(m)
print(1/m)
#may need to remove ends, should use jupyter notebook or store the data somewhereo

#As there is no preference to different points on the strange attractor I think that this is viable


"""stuff to do:
-try and get liapunov coefficent from sensativity calculation
-try and calculate the dimention
-vary r and look at the r=24.78
-also look at case beforehand (phase protraits may be nice here)
-write up clearly what you find and comment code making it clear what you have done (may be useful to
store data in jsons so you don't have to keep running it)
 
 """