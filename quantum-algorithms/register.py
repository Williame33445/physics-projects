import numpy as np
import copy


q0 = np.array([1, 0],dtype=np.complex128)  # |0>
q1 = np.array([0, 1],dtype=np.complex128)  # |1>

def intialiseProductRegister(productOfQubits):
    qubits = np.array([1],dtype=np.complex128)
    for q in productOfQubits:
        qNorm = q/np.linalg.norm(q)
        qubits = np.kron(qubits, qNorm)
    return Register(qubits)

def mergeRegisters(r1,r2):
    return Register(np.kron(r1.qubits,r2.qubits))

def mergeRegisterOperator(o1,o2):
    return np.kron(o1,o2)

class Register:
    def __init__(self, qubits): #this needs to be altered to a more general situation
        self.N = np.shape(qubits)[0]
        self.n = np.log2(self.N).astype(int)

        self.qubits = qubits
        
        self.binary = []
        for i in range(self.N):
            self.binary.append(format(i, f'0{self.n}b'))


    def displayQubits(self,dType="int",cutoff=3,binBits=0,plusNumber=20):
        registerString = ""
        for i in range(self.N):
            if np.round(self.qubits[i],cutoff) != 0:
                if dType == "int":
                    index = str(i)
                elif dType == "bin":
                    index = self.binary[i]
                elif dType == "intbin":
                    index = str(int(self.binary[i][:-binBits],2)) + "," + self.binary[i][-binBits:]
                elif dType == "binint":
                    index = self.binary[i][:-binBits] + "," + str(int(self.binary[i][-binBits:],2))

                registerString += f" + {np.round(self.qubits[i],cutoff)}|" + index + ">"
        #only return the first plusNumber + terms
        splitString = registerString.split("+")
        if len(splitString) <= plusNumber:  
            return registerString[3:]
        return "+".join(registerString[3:].split("+")[:plusNumber] + [" ..."])
    
    def findSingleQubitOp(self,op,index):
        registerOp = np.array([1],dtype=np.complex128)
        for i in range(self.n): 
            if i == index:
                registerOp = np.kron(registerOp, op)
            else:
                registerOp = np.kron(registerOp, np.eye(2,dtype=np.complex128))
        return registerOp

    def findControlledSingleQubitOp(self,op, controlIndex, targetIndex):
        
        nonControlled = self.findSingleQubitOp(np.array([[1,0],[0,0]],dtype=np.complex128), controlIndex)
        controlled = self.findSingleQubitOp(np.array([[0,0],[0,1]],dtype=np.complex128), controlIndex) @ self.findSingleQubitOp(op, targetIndex)
        return nonControlled + controlled
    
    def getSingleQubitProduct(self,op,start,stop):
        registerOp = np.array([1],dtype=np.complex128)
        for i in range(self.n):
            if start <= i < stop:
                registerOp = np.kron(registerOp, op)
            else:
                registerOp = np.kron(registerOp, np.eye(2,dtype=np.complex128))
        return registerOp

    def applySingleQubitProduct(self,op,start,stop):
        self.qubits = self.getSingleQubitProduct(op,start,stop) @ self.qubits

        
    def applyGate(self, gate):
        self.qubits = gate @ self.qubits

    def applyH(self, index):
        H = np.array([[1, 1], [1, -1]],dtype=np.complex128) / np.sqrt(2,dtype=np.complex128)
        self.qubits = self.findSingleQubitOp(H, index) @ self.qubits

    def applyCNOT(self, controlIndex, targetIndex):
        X = np.array([[0,1],[1,0]],dtype=np.complex128)
        CNOT = self.findControlledSingleQubitOp(self, X, controlIndex, targetIndex)
        self.qubits = CNOT @ self.qubits

    def invertQubits(self,start,stop):
        reveresedBinary = [int(b[:start] + b[start:stop][::-1] + b[stop:],2) for b in self.binary]
        newRegister = copy.deepcopy(self.qubits)
        for i in range(len(reveresedBinary)):
            newRegister[reveresedBinary[i]] = self.qubits[i]
        self.qubits = newRegister     