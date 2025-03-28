import numpy as np

q0 = np.array([1, 0],dtype=np.complex_)  # |0>
q1 = np.array([0, 1],dtype=np.complex_)  # |1>

class Register:
    def __init__(self, qubits, cutoff=10**-9,numberType="int"):
        self.n = len(qubits)
        self.N = 2**self.n

        self.cutoff = cutoff
        self.numberType = numberType

        self.register = np.array([1],dtype=np.complex_) #better name given what the class is called
        for q in qubits:
            qNorm = q/np.linalg.norm(q)
            self.register = np.kron(self.register, qNorm)

    def __str__(self):
        registerString = ""
        for i in range(self.N):
            if np.linalg.norm(self.register[i]) > self.cutoff:
                if self.numberType == "int":
                    index = str(i)
                elif self.numberType == "bin":
                    index = format(i, f'0{self.n}b')

                registerString += f" + {np.round(self.register[i],3)}|" + index + ">"
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
    
    def findSingleQubitProduct(self,op):
        registerOp = np.array([1],dtype=np.complex_)
        for i in range(self.n):
            registerOp = np.kron(registerOp, op)
        return registerOp
    
    def applyGate(self, gate):
        self.register = gate @ self.register

    def applyH(self, index): #this combined with TOFFOLI gate gives universal quantum computation
        H = np.array([[1, 1], [1, -1]],dtype=np.complex_) / np.sqrt(2,dtype=np.complex_)
        op = self.findSingleQubitOp(H, index)
        self.register = op @ self.register

    def applyCNOT(self, controlIndex, targetIndex):
        X = np.array([[0,1],[1,0]],dtype=np.complex_)
        CNOT = self.findControlledSingleQubitOp(self, X, controlIndex, targetIndex)
        self.register = CNOT @ self.register

    def applyTOFFOLI(self, controlIndex1, controlIndex2, targetIndex):
        #need to think about this 
        pass
        

        
# add measurement part as well

#QFT example
QFTRegister = Register([q0,q0,q1])
print(f"Initial state: {str(QFTRegister)}")

def R_m(m):
    return np.array([[1,0],[0,np.exp(2*np.pi*1j/(2**m),dtype=np.complex_)]],dtype=np.complex_)

Rs = [R_m(i) for i in range(2,QFTRegister.n+1)]
for i in range(QFTRegister.n):
    QFTRegister.applyH(i)

    for j in range(QFTRegister.n - i - 1):
        gate = QFTRegister.findControlledSingleQubitOp(Rs[j], i+1+j, i)
        QFTRegister.applyGate(gate)

print(f"Final state: {str(QFTRegister)}")