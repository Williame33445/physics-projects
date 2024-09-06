import basis_set_exchange as basis_set_exchange 
import numpy as np

#code to get GTO data from Basis Stack exchange

class GTOPrimative:
    """
    Class that holds the primative GTO.
    """
    def __init__(self,alpha,center,type):
        self.alpha = alpha
        self.center = center
        self.type = type
    
class GTOContraction:
    """
    Class that holds the contracted GTOs.
    """
    def __init__(self,primitives,ds):
        self.primitives = primitives 
        self.ds = np.array(ds)
        
def getType(l):
    """
    Finds all GTO types for a given l.
    """
    if l == 0:
        return [np.array([0,0,0])]
    elif l == 1:
        return [np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])]
    elif l == 2:
        return [np.array([1,1,0]),np.array([1,0,1]),np.array([0,1,1]),np.array([2,0,0]),np.array([0,2,0])]
    else:
        print("Error, l required is greater than implemented in getType")

def getGuassians(set,atom,center,basisType):
    """
    Gets a certain Gaussian set from Basis Stack exchange and puts it into a "contracted" or "primative"
    form.
    """
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
