import basis_set_exchange as basis_set_exchange 
import numpy as np
from GTOs import GTOPrimitive, GTOContraction

#code to get GTO data from Basis Stack exchange

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
    elif l == 3:
        return [np.array([3,0,0]),np.array([0,3,0]),np.array([0,1,2]),np.array([2,1,0]),np.array([0,2,1]),
                np.array([1,2,0]),np.array([1,0,2])]
    else:
        print("Error, l required is greater than implemented in getType")

def getGuassians(set,atom,center,basisType):
    """
    Gets a certain Gaussian set from Basis Stack exchange and puts it into a "contracted" or "primitive"
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
                        contractedGTOs.append(GTOPrimitive(float(exponents[i]),center,type))
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
