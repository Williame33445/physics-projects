from scipy.linalg import eigh

class Electron:
    """
    Class that holds the electron state.
    """
    def __init__(self,energy,ds,spin):
        self.energy = energy
        self.ds = ds
        self.spin = spin


def takeGroundEigStates(sortedStates,N):
    """
    Function that returns the states corresponding to the ground state for N electrons.
    """
    return sortedStates[:N]

def takeLthExcitedSpin(sortedStates,N,L,spin):
    """
    Function that returns the electrons corresponding to a certain spin being excited to its Lth
    state (N is the number of electrons).
    """
    spinStates = list(filter(lambda s: s.spin == spin, sortedStates[N-1:]))
    return sortedStates[:N-1] + [spinStates[L]]

     

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

    #put eigenstates into electron objects and combine
    combinedUp = [Electron(energy,eigenVecsUp[:,i],"up") for i,energy in enumerate(eigenValsUp)]
    combinedDown = [Electron(energy,eigenVecsDown[:,i],"down") for i,energy in enumerate(eigenValsDown)]
    combinedStates = combinedUp + combinedDown

    #find target eigenstates (usually ground) and normalise
    sortedStates = sorted(combinedStates,key=lambda s: s.energy)
    ocuppiedStates = rep.normaliseElectronList(getTargetEigStates(sortedStates)) 

    #find energy
    ENew = rep.findE(ocuppiedStates)

    #is result is accurate enough/ recursion depth is reached return, otherwise iterate again
    if abs(ENew - E) < maxError or depth>150:
        return ENew,ocuppiedStates
    elif depth>999:
        print("max recursion depth")
        return ENew,ocuppiedStates
    else:
        upsNew = [s.ds for s in ocuppiedStates if s.spin == "up"]
        downsNew = [s.ds for s in ocuppiedStates if s.spin == "down"]
        return iterateHF(upsNew,downsNew,rep,ENew,maxError,getTargetEigStates,depth+1)
