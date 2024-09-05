import basis_set_exchange as basis_set_exchange 
import numpy as np

#code to get GTO data from Basis Stack exchange

class GTOPrimative:
    def __init__(self,alpha,center,type):
        self.alpha = alpha
        self.center = center
        self.type = type

class GTOContraction:
    def __init__(self,primitives,ds):
        self.primitives = primitives 
        self.ds = np.array(ds)
        
def getType(l):
    if l == 0:
        return [np.array([0,0,0])]
    elif l == 1:
        return [np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])]
    elif l == 2:
        return [np.array([1,1,0]),np.array([1,0,1]),np.array([0,1,1]),np.array([2,0,0]),np.array([0,2,0])]
    else:
        print("over")

def getGuassians(set,atom,center,basisType):
    bs_dict = basis_set_exchange.get_basis(set,elements=atom)
    GTOs = []
    for e in bs_dict["elements"].values():
        for contraction in e["electron_shells"]:
            contractedGTOs = []
            ds = []
            for l in contraction["angular_momentum"]:
                typeL = getType(l)
                for type in typeL:
                    exponents = contraction["exponents"]
                    coefficients = contraction["coefficients"]
                    for i in range(len(exponents)):
                        contractedGTOs.append(GTOPrimative(float(exponents[i]),center,type))
                        #this is to take into account Pople type basis sets, may cause error in some cases
                        try: 
                            ds.append(float(coefficients[l][i]))
                        except:
                            #print("Assumed Pople type conventions are followed")
                            ds.append(float(coefficients[0][i]))
            
            if basisType == "primitive":
                GTOs += contractedGTOs
            elif basisType == "contracted":
                GTOs.append(GTOContraction(contractedGTOs,ds))
            else:
                print("Type not formated correctly")
                
    return GTOs
