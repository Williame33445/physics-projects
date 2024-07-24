import numpy as np

#data for potential
h = 0.005
rs = np.arange(0,0.5,step=h)
l = 8

def f(r):
    if r == 0:
        return 0
    return (l*(l+1)/r**2)

us = np.empty(len(rs))
us[0] = 0
us[1] = h**(l+1)

def w(u,r):
    return (1 - f(r)*(h**2)/12)*u
def u(w,r):
    return w/(1 - f(r)*(h**2)/12)

ws = np.vectorize(w)(us,rs)
for n in range(1,len(rs)-1):
    ws[n+1] = 2*ws[n] - ws[n-1] + (h**2)*us[n]*f(rs[n])
    us[n+1] = u(ws[n+1],rs[n+1])

    

t1 = rs**(l+1)
t2 = us#*t1[-1]/us[-1] 
print(t2[3])
print(t1[3])

print(abs(t1-t2))

