import numpy as np

class Eigenstates:
    def __init__(self,ups,downs,rep):
        self.ups = ups
        self.downs = downs
        self.rep = rep

        for i in range(len(self.ups)):
            self.ups[i] = self.rep.normalise(self.ups[i])
            self.downs[i] = self.rep.normalise(self.downs[i])

    def getCs(self,part): #is there a nicer way to do this?
        if part == "up":
            return self.ups
        elif part == "down":
            return self.downs

    def getOppositeCs(self,part):
        if part == "up":
            return self.downs
        elif part == "down":
            return self.ups

    def F(self,spinDirection):
        F = np.empty([self.rep.repNumber,self.rep.repNumber])
        for p in range(self.rep.repNumber):
            for q in range(self.rep.repNumber):
                F[p,q] = self.rep.h[p,q]

                for CList in self.getCs(spinDirection):
                    for r in range(self.rep.repNumber): #can be made more efficient
                        for s in range(self.rep.repNumber):
                            F[p,q] += CList[r]*CList[s]*(self.rep.twoElecInts[p,r,q,s] - self.rep.twoElecInts[p,r,s,q])

                
                for CList in self.getOppositeCs(spinDirection):
                    for r in range(self.rep.repNumber):
                        for s in range(self.rep.repNumber): 
                            F[p,q] += CList[r]*CList[s]*self.rep.twoElecInts[p,r,q,s]
        
        return F
    
    def findE(self,epsilons):
        E = sum(epsilons)
        for C in self.ups:
            for i,c_i in enumerate(C):
                for j,c_j in enumerate(C): 
                    E += self.rep.h[i,j]*c_i*c_j
        
        for C in self.downs:
            for i,c_i in enumerate(C):
                for j,c_j in enumerate(C):
                    E += self.rep.h[i,j]*c_i*c_j
        
        return E/2
    
    def solveEigenEquation(self,spinDirection):
        return np.linalg.eig(np.matmul(self.rep.invS,self.F(spinDirection)))
        
