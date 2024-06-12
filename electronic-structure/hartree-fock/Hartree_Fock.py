import numpy as np
#why do C's have to be real?

class Eigenstates:
    def __init__(self,ups,downs):
        self.ups = ups
        self.downs = downs

    def getCs(self,part):
        if part == "up":
            return self.ups
        elif part == "down":
            return self.downs
    def getOppositeCs(self,part):
        if part == "up":
            return self.downs
        elif part == "down":
            return self.ups


class HF:
    def __init__(self,Zs,N,basis,eigenstatesGuess,EGuess,maxError):
        self.Zs = Zs
        self.N = N #number of electrons
        self.basis = basis #is this the correct name?
        self.eigenstates = eigenstatesGuess # might be better to normalise somewhere else
        self.E = EGuess
        self.maxError = maxError

        self.eigenstates.ups = [self.normalise(i) for i in self.eigenstates.ups]
        self.eigenstates.downs = [self.normalise(i) for i in self.eigenstates.downs]


        self.iterateHF()

    def solveEigenEquation(self,spinDirection):
        return np.linalg.eig(np.matmul(np.linalg.inv(self.basis.S),self.F(spinDirection)))
    
    def findE(self,epsilons):
        E = sum(epsilons)
        for C in self.eigenstates.ups:
            for i,c_i in enumerate(C):
                for j,c_j in enumerate(C): 
                    E += self.basis.h[i,j]*c_i*c_j
        
        for C in self.eigenstates.downs:
            for i,c_i in enumerate(C):
                for j,c_j in enumerate(C):
                    E += self.basis.h[i,j]*c_i*c_j
        
        return E/2


        
    def iterateHF(self):
        eigenvaluesUp,eigenstatesUp = self.solveEigenEquation("up")
        eigenvaluesDown,eigenstatesDown = self.solveEigenEquation("down")

        #is there a better way to do this bit?

        combineParamsUp = lambda eigenvalue,eigenstate: {"epsilon": eigenvalue, "state": self.normalise(eigenstate), "spin":"up"}
        combinedUp = list(map(combineParamsUp,eigenvaluesUp,eigenstatesUp))

        combineParamsDown = lambda eigenvalue,eigenstate: {"epsilon": eigenvalue, "state": self.normalise(eigenstate), "spin":"down"}
        combinedDown = list(map(combineParamsDown,eigenvaluesDown,eigenstatesDown))

        combinedStates = combinedUp + combinedDown

        sortedStates = sorted(combinedStates,key=lambda v: v["epsilon"])
        groundStates = sortedStates[:self.N]

        self.eigenstates.ups = [s["state"] for s in groundStates if s["spin"] == "up"]
        self.eigenstates.downs = [s["state"] for s in groundStates if s["spin"] == "down"]
        

        E = self.findE([s["epsilon"] for s in groundStates])

        if abs(E - self.E) < self.maxError:
            print(E)
            pass
        else:
            self.E = E
            self.iterateHF()










        
    def normalise(self,C):
        innerProduct = 0
        for i,c_i in enumerate(C):
            for j,c_j in enumerate(C):
                innerProduct += c_i*c_j*self.basis.S[i][j]
        
        
        return C/(innerProduct**0.5)

        

    def F(self,spinDirection):
        F = np.empty([self.basis.basisNumber,self.basis.basisNumber])
        for p in range(self.basis.basisNumber):
            for q in range(self.basis.basisNumber):
                F[p,q] = self.basis.h[p,q]

                for CList in self.eigenstates.getCs(spinDirection):
                    for r in range(self.basis.basisNumber): #can be made more efficient
                        for s in range(self.basis.basisNumber):
                            F[p,q] += CList[r]*CList[s]*(self.basis.twoElecInts[p,r,q,s] - self.basis.twoElecInts[p,r,s,q])

                
                for CList in self.eigenstates.getOppositeCs(spinDirection):
                    for r in range(self.basis.basisNumber):
                        for s in range(self.basis.basisNumber): 
                            F[p,q] += CList[r]*CList[s]*self.basis.twoElecInts[p,r,q,s]
        
        return F





