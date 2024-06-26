#method solves d^2X/dt^2 = F(X,t) to h^2

import numpy as np

def iterateVerlet(x_tPlush,x_t,h,tPlush,F):
    #F = F(x,t)
    return 2*x_tPlush - x_t + F(x_tPlush,tPlush)*h**2
