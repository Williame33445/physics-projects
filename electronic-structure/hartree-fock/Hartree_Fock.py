import numpy as np
#why do C's have to be real?

#this doesn't need to be a class

def F(eigenstates1,eigenstates2,rep): #better names
    F = np.empty([rep.repNumber,rep.repNumber])
    for p in range(rep.repNumber):
        for q in range(rep.repNumber):
            F[p,q] = rep.h[p,q]

            for C in eigenstates1:
                for r,c_r in enumerate(C): #can be made more efficient
                    for s,c_s in enumerate(C):
                        F[p,q] += c_r*c_s*(rep.twoElecInts[p,r,q,s] - rep.twoElecInts[p,r,s,q])

            
            for C in eigenstates2:
                for r,c_r in enumerate(C):
                    for s,c_s in enumerate(C): 
                        F[p,q] += c_r*c_s*rep.twoElecInts[p,r,q,s]
    
    return F
    
    
        


def findE(occupiedStates,rep):
    E = 0
    for S in occupiedStates:
        E += S["epsilon"]
        for i,c_i in enumerate(S["state"]):
            for j,c_j in enumerate(S["state"]):
                E += rep.h[i,j]*c_i*c_j
    
    return E/2



def iterateHF(ups,downs,rep,N,E,maxError):
    FPlus = F(ups,downs,rep) #ups and downs need to be normalised before
    FMinus = F(downs,ups,rep)


    eigenvaluesUp,eigenstatesUp = np.linalg.eig(np.matmul(rep.invS,FPlus))
    eigenvaluesDown,eigenstatesDown = np.linalg.eig(np.matmul(rep.invS,FMinus))

    combineParamsUp = lambda eigenvalue,eigenstate: {"epsilon": eigenvalue, "state": eigenstate, "spin":"up"}
    combinedUp = list(map(combineParamsUp,eigenvaluesUp,eigenstatesUp))

    combineParamsDown = lambda eigenvalue,eigenstate: {"epsilon": eigenvalue, "state": eigenstate, "spin":"down"}
    combinedDown = list(map(combineParamsDown,eigenvaluesDown,eigenstatesDown))

    combinedStates = combinedUp + combinedDown

    sortedStates = sorted(combinedStates,key=lambda v: v["epsilon"])
    ocuppiedStates = rep.normaliseDicList(sortedStates[:N]) #generalise

    upsNew = [s["state"] for s in ocuppiedStates if s["spin"] == "up"]
    downsNew = [s["state"] for s in ocuppiedStates if s["spin"] == "down"]

    
    ENew = findE(ocuppiedStates,rep)
    

    if abs(ENew - E) < maxError:
            print(ENew)
    else:
        upsNew = [s["state"] for s in ocuppiedStates if s["spin"] == "up"]
        downsNew = [s["state"] for s in ocuppiedStates if s["spin"] == "down"]
        iterateHF(upsNew,downsNew,rep,N,ENew,maxError)