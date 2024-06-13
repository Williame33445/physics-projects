import numpy as np
#why do C's have to be real?

def combineParams(eigenVal,eigenVec,spin):
    return {"e": eigenVal, "state": eigenVec, "spin": spin}


def iterateHF(ups,downs,rep,N,E,maxError):
    FPlus = rep.F(ups,downs) #ups and downs need to be normalised before
    FMinus = rep.F(downs,ups)


    eigenValUp,eigenVecUp = np.linalg.eig(np.matmul(rep.invS,FPlus))
    eigenValDown,eigenVecDown = np.linalg.eig(np.matmul(rep.invS,FMinus))

    combinedUp = list(map(lambda val,vec:{"e": val, "state": vec, "spin":"up"},eigenValUp,eigenVecUp))
    combinedDown = list(map(lambda val,vec:{"e": val, "state": vec, "spin":"down"},eigenValDown,eigenVecDown))
    combinedStates = combinedUp + combinedDown

    sortedStates = sorted(combinedStates,key=lambda v: v["e"])
    ocuppiedStates = rep.normaliseDicList(sortedStates[:N]) #generalise
    
    ENew = rep.findE(ocuppiedStates)
    

    if abs(ENew - E) < maxError:
            return ENew
    else:
        upsNew = [s["state"] for s in ocuppiedStates if s["spin"] == "up"]
        downsNew = [s["state"] for s in ocuppiedStates if s["spin"] == "down"]
        return iterateHF(upsNew,downsNew,rep,N,ENew,maxError)
