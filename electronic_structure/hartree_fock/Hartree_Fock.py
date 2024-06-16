import numpy as np
from scipy.linalg import eigh

#why do C's have to be real?

def takeGroundEigStates(sortedStates,N):
    return sortedStates[:N]
     

def iterateHF(ups,downs,rep,E,maxError,getTargetEigStates,depth=0): 
    FPlus = rep.F(ups,downs) #ups and downs need to be normalised before
    FMinus = rep.F(downs,ups)


    eigenValsUp,eigenVecsUp = eigh(FPlus,rep.S)
    eigenValsDown,eigenVecsDown = eigh(FMinus,rep.S)

    combinedUp = [{"e":val,"state":eigenVecsUp[:,i],"spin":"up"} for i,val in enumerate(eigenValsUp)]
    combinedDown = [{"e":val,"state":eigenVecsDown[:,i],"spin":"down"} for i,val in enumerate(eigenValsDown)]
    combinedStates = combinedUp + combinedDown

    sortedStates = sorted(combinedStates,key=lambda v: v["e"])
    ocuppiedStates = rep.normaliseDicList(getTargetEigStates(sortedStates)) 
    
    ENew = rep.findE(ocuppiedStates)
    

    if abs(ENew - E) < maxError or depth>150:
        return ENew,ocuppiedStates
    elif depth>600:
        print("max recursion depth")
        return ENew,ocuppiedStates
    else:
        upsNew = [s["state"] for s in ocuppiedStates if s["spin"] == "up"]
        downsNew = [s["state"] for s in ocuppiedStates if s["spin"] == "down"]
        return iterateHF(upsNew,downsNew,rep,ENew,maxError,getTargetEigStates,depth+1)
