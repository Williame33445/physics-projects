from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import numpy as np
from spherically_symetric_quantum_scattering import ScatteringSystem

#This program uses the ScatteringSystem class to calculate the distribution of total cross section 
# vs energy for the K-H scattering.


#constants, energy is measured in meV and distance in units of rho
epsilon = 5.9
rho = 1
C = np.sqrt(epsilon*6.12/25) #needed for initial conditions

def V_LJ(r):
    """
    Lennard-Jones potential.
    """
    return epsilon*((rho/r)**12 - 2*((rho/r)**6))

#this function is needed as 2 intial boundary conditions are required for the numerical simulation 
def uSmallr(r):
    """
    Approximate solution to u for small r
    """
    return np.exp(-C/(r**5))


#find intial conditions and h 
r_0 = 0.7*rho
r_1 = 0.701*rho
h = r_1 - r_0
u_0 = uSmallr(r_0)
u_1 = uSmallr(r_1)

#system parameters
EArray = np.arange(.1,3.5,.01)
r_end = 5

#run simulation
totalCrossSectionList = []
for i,E in enumerate(EArray):
    totalCrossSectionList.append(ScatteringSystem(E,V_LJ,r_0,h,u_0,u_1,r_end,8).totalCrossSection)
    print(str(np.round(100*i/len(EArray),decimals=2))+"%")

#graph
plt.plot(EArray,totalCrossSectionList)
plt.xlabel("Energy [MeV]")
plt.ylabel(r'Total Cross Section [$\rho^2$]')
plt.title("Total Cross section vs energy for H-Kr scattering")
plt.show()