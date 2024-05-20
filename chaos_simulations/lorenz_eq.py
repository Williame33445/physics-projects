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

#simulate function
def simulate(RStates):
    t = 0
    for x in range(10**5):
        RStates.append(iterateRK4(t,RStates[-1],h,f))
    return [[s[c] for s in RStates] for c in range(3)]


# example simulation and graph
states = simulate([np.array([0.1,0,0])])

ax = plt.axes(projection='3d')
ax.plot3D(states[0], states[1], states[2], 'blue',linewidth = '.5')
plt.show()

# Use this to generate a large number of simulations (via GPUs?) then try and get dim via averages
num = 2
evenList = []
for x in np.linspace(-.5,.5,num):
    for y in np.linspace(-.5,.5,num):
        for z in np.linspace(-.5,.5,num):
            evenList.append(np.array([x,y,z]))

#As there is no preference to different points on the strange attractor I think that this is viable


"""stuff to do:
-look at the sensitivity on intial conditions (deviation of close solutions), try and find exponential
divergence
-try and calculate the dimention
-vary r and look at the r=24.78
-also look at case beforehand (phase protraits may be nice here)
-write up clearly what you find and comment code making it clear what you have done (may be useful to
store data in jsons so you don't have to keep running it)
 
 """