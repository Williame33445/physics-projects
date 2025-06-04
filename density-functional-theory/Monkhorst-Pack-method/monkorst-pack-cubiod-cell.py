#let lattice constants be a1,a2,a3
#let sampling range be q1,q2,q3 (have to be natural numbers)

import numpy as np

def u(r,q):
    return (2*r - q - 1)/(2*q)

def MPIntegrateCubiod(f, q1, q2, q3, a1, a2, a3):
    b1 = 2*np.pi*np.array([1/a1,0,0])
    b2 = 2*np.pi*np.array([0,1/a2,0])
    b3 = 2*np.pi*np.array([0,0,1/a3])
    
    integral = 0
    for r1 in range(1,q1+1):
        for r2 in range(1,q2+1):
            for r3 in range(1,q3+1):
                integral += f(u(r1,q1)*b1 + u(r2,q2)*b2 + u(r3,q3)*b3)/(q1*q2*q3)
    return integral


#test, exactly the integral is ((pi/a)^2)/3
def fTest(k):
    return k[0]**2

print("If BZ samples in each direction is choosen as 2: ")
for a in range(1,5):
    print(f"Size: {a}. Error: {(MPIntegrateCubiod(fTest, 2, 2, 2, a, a, a)-((np.pi/a)**2/3))}")

#This is what we expect as the size is increased the absolute error decreases, as we are using the 
#same number of samples for a smaller volume. (note that the relative error remains the same as expected)

print("If size of the system is choosen as 1: ")
for s in range(1,5):
    print(f"Sample size in each direction: {s}. Error: {(MPIntegrateCubiod(fTest, s, s, s, 1, 1, 1)-(np.pi**2/3))}")

#As expected accuracy increases with the number of samples