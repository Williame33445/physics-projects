import sys
import os

sys.path.append(os.path.abspath("."))
from useful_code.RK4_methods import *
from mpl_toolkits.mplot3d import Axes3D


import numpy as np
import matplotlib.pyplot as plt


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

#actual simulation
h = 0.001
t = 0
RStates = [np.array([0.1,0,0])]
for x in range(10**5):
    RStates.append(iterateRK4(t,RStates[-1],h,f))


#graph
ax = plt.axes(projection='3d')

xline = [i[0] for i in RStates]
yline = [i[1] for i in RStates]
zline = [i[2] for i in RStates]
ax.plot3D(xline, yline, zline, 'blue',linewidth = '.5')
plt.show()


