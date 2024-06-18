import numpy as np
from scipy.linalg import eigh


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
    combinedUp = [{"e":val,"state":eigenVecsUp[:,i],"spin":"up"} for i,val in enumerate(eigenValsUp)]
    combinedDown = [{"e":val,"state":eigenVecsDown[:,i],"spin":"down"} for i,val in enumerate(eigenValsDown)]
    combinedStates = combinedUp + combinedDown

    #find target eigenstates (usually ground) and normalise
    sortedStates = sorted(combinedStates,key=lambda v: v["e"])
    ocuppiedStates = rep.normaliseDicList(getTargetEigStates(sortedStates)) 

    #find energy
    ENew = rep.findE(ocuppiedStates)

    #is result is accurate enough/ recursion depth is reached return, otherwise iterate again
    if abs(ENew - E) < maxError or depth>150:
        return ENew,ocuppiedStates
    elif depth>600:
        print("max recursion depth")
        return ENew,ocuppiedStates
    else:
        upsNew = [s["state"] for s in ocuppiedStates if s["spin"] == "up"]
        downsNew = [s["state"] for s in ocuppiedStates if s["spin"] == "down"]
        return iterateHF(upsNew,downsNew,rep,ENew,maxError,getTargetEigStates,depth+1)
