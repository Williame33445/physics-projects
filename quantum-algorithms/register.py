import numpy as np
import copy


q0 = np.array([1, 0],dtype=np.complex_)  # |0>
q1 = np.array([0, 1],dtype=np.complex_)  # |1>

class Register:
    def __init__(self, productOfQubits):
        self.n = len(productOfQubits)
        self.N = 2**self.n
        
        self.qubits = np.array([1],dtype=np.complex_)
        for q in productOfQubits:
            qNorm = q/np.linalg.norm(q)
            self.qubits = np.kron(self.qubits, qNorm)

        self.binary = []
        for i in range(self.N):
            self.binary.append(format(i, f'0{self.n}b'))


    def displayQubits(self,cutoff=10**-9,endBitsRemoved=0,type="int"):
        registerString = ""
        for i in range(self.N):
            if np.linalg.norm(self.qubits[i]) > cutoff:
                if endBitsRemoved > 0:
                    index = self.binary[i][:-endBitsRemoved]
                else:
                    index = self.binary[i]

                if type == "int":
                    index = str(int(index,2))

                registerString += f" + {np.round(self.qubits[i],3)}|" + index + ">"
        return registerString[3:]
    
    def findSingleQubitOp(self,op,index):
        registerOp = np.array([1],dtype=np.complex_)
        for i in range(self.n): 
            if i == index:
                registerOp = np.kron(registerOp, op)
            else:
                registerOp = np.kron(registerOp, np.eye(2,dtype=np.complex_))
        return registerOp

    def findControlledSingleQubitOp(self,op, controlIndex, targetIndex):
        
        nonControlled = self.findSingleQubitOp(np.array([[1,0],[0,0]],dtype=np.complex_), controlIndex)
        controlled = self.findSingleQubitOp(np.array([[0,0],[0,1]],dtype=np.complex_), controlIndex) @ self.findSingleQubitOp(op, targetIndex)
        return nonControlled + controlled
    
    def applySingleQubitProduct(self,op,start,stop):
        registerOp = np.array([1],dtype=np.complex_)
        for i in range(self.n):
            if start <= i < stop:
                registerOp = np.kron(registerOp, op)
            else:
                registerOp = np.kron(registerOp, np.eye(2,dtype=np.complex_))
        self.qubits = registerOp @ self.qubits
        
    def applyGate(self, gate):
        self.qubits = gate @ self.qubits

    def applyH(self, index):
        H = np.array([[1, 1], [1, -1]],dtype=np.complex_) / np.sqrt(2,dtype=np.complex_)
        self.qubits = self.findSingleQubitOp(H, index) @ self.qubits

    def applyCNOT(self, controlIndex, targetIndex):
        X = np.array([[0,1],[1,0]],dtype=np.complex_)
        CNOT = self.findControlledSingleQubitOp(self, X, controlIndex, targetIndex)
        self.qubits = CNOT @ self.qubits

    def invertQubits(self,start,stop):
        reveresedBinary = [int(b[:start] + b[start:stop][::-1] + b[stop:],2) for b in self.binary]
        newRegister = copy.deepcopy(self.qubits)
        for i in range(len(reveresedBinary)):
            newRegister[reveresedBinary[i]] = self.qubits[i]
        self.qubits = newRegister     
