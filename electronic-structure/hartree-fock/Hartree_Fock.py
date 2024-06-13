import numpy as np
from eigenstates import *
#why do C's have to be real?

#this doesn't need to be a class
class HF:
    def __init__(self,Zs,N,rep,eigenstatesGuess,EGuess,maxError):
        self.Zs = Zs
        self.N = N #number of electrons
        self.rep = rep #is this the correct name?
        self.eigenstates = eigenstatesGuess # might be better to normalise somewhere else
        self.E = EGuess
        self.maxError = maxError


        self.iterateHF()

    


        
    def iterateHF(self):
        eigenvaluesUp,eigenstatesUp = self.eigenstates.solveEigenEquation("up")
        eigenvaluesDown,eigenstatesDown = self.eigenstates.solveEigenEquation("down") #might want to do this automatically and create new eigenstate each time?

        #is there a better way to do this bit?

        combineParamsUp = lambda eigenvalue,eigenstate: {"epsilon": eigenvalue, "state": eigenstate, "spin":"up"}
        combinedUp = list(map(combineParamsUp,eigenvaluesUp,eigenstatesUp))

        combineParamsDown = lambda eigenvalue,eigenstate: {"epsilon": eigenvalue, "state": eigenstate, "spin":"down"}
        combinedDown = list(map(combineParamsDown,eigenvaluesDown,eigenstatesDown))

        combinedStates = combinedUp + combinedDown

        sortedStates = sorted(combinedStates,key=lambda v: v["epsilon"])
        groundStates = sortedStates[:self.N]

        ups = [s["state"] for s in groundStates if s["spin"] == "up"]
        downs = [s["state"] for s in groundStates if s["spin"] == "down"]

        self.eigenstates = Eigenstates(ups,downs,self.rep)
        

        E = self.eigenstates.findE([s["epsilon"] for s in groundStates])

        if abs(E - self.E) < self.maxError:
            print(E)
            pass
        else:
            self.E = E
            self.iterateHF()
