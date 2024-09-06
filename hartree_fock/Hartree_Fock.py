import numpy as np
from scipy.linalg import eigh

class Shell:
    def __init__(self,energy,ds,spin):
        self.energy = energy
        self.ds = ds
        self.spin = spin


def takeGroundEigStates(sortedStates,N):
    """
    Function that returns the states corresponding to the ground state for N electrons.
    """
    return sortedStates[:N]
     

def iterateHF(ups,downs,rep,E,maxError,getTargetEigStates,depth=0): 
    """
    Applies linear Hatree-Fock iteration.

    Inputs:
    ups - C^+_k matrix from (5)
    downs - C^-_k matrix from (5)
    rep - representation class (contians matrix elements)
    E - Energy found in previous iteration
    getTargetEigStates - function that returns target eigenstate (eg. ground state)
    depth - int that measures the recursion depth
    """
    #Find F s
    FPlus = rep.F(ups,downs)
    FMinus = rep.F(downs,ups)

    #find eigenvalues, eigh solves generalised eigenequation
    eigenValsUp,eigenVecsUp = eigh(FPlus,rep.S)
    eigenValsDown,eigenVecsDown = eigh(FMinus,rep.S)

    #label eigenstates as dictionaries and combine
    combinedUp = [Shell(energy,eigenVecsUp[:,i],"up") for i,energy in enumerate(eigenValsUp)]
    combinedDown = [Shell(energy,eigenVecsDown[:,i],"down") for i,energy in enumerate(eigenValsDown)]
    combinedStates = combinedUp + combinedDown

    #find target eigenstates (usually ground) and normalise
    sortedStates = sorted(combinedStates,key=lambda s: s.energy)
    ocuppiedStates = rep.normaliseShell(getTargetEigStates(sortedStates)) 

    #find energy
    ENew = rep.findE(ocuppiedStates)

    #is result is accurate enough/ recursion depth is reached return, otherwise iterate again
    if abs(ENew - E) < maxError or depth>150:
        return ENew,ocuppiedStates
    elif depth>600:
        print("max recursion depth")
        return ENew,ocuppiedStates
    else:
        upsNew = [s.ds for s in ocuppiedStates if s.spin == "up"]
        downsNew = [s.ds for s in ocuppiedStates if s.spin == "down"]
        return iterateHF(upsNew,downsNew,rep,ENew,maxError,getTargetEigStates,depth+1)
