import numpy as np

class GTOPrimitive:
    """
    Class that holds the primitive GTO.
    """
    def __init__(self,alpha,center,type):
        self.alpha = alpha
        self.center = center
        self.type = type
    
    def wavefunction(self,xv,yv,zv):
        Rx = self.center[0]
        Ry = self.center[1]
        Rz = self.center[2]
        i = self.type[0]
        j = self.type[1]
        k = self.type[2]
        return (xv - Rx)**i*(yv - Ry)**j*(zv - Rz)**k*np.exp(-self.alpha*((xv - Rx)**2 + (yv - Ry)**2 + (zv - Rz)**2))
    
    
class GTOContraction:
    """
    Class that holds the contracted GTOs.
    """
    def __init__(self,primitives,ds):
        self.primitives = primitives 
        self.ds = ds
    
    def wavefunction(self,xv,yv,zv):
        w = np.zeros(xv.shape)
        for i,d in enumerate(self.ds):
            w += d*self.primitives[i].wavefunction(xv,yv,zv)
        return w
        
def getDensity(xv,yv,zv,electrons,GTOs):
    GTOWaveFunc = []
    for GTO in GTOs:
        GTOWaveFunc.append(GTO.wavefunction(xv,yv,zv))

    ns = np.zeros(xv.shape)
    for electron in electrons:
        ne = np.zeros(xv.shape)
        for i,d in enumerate(electron.ds):
            ne += d*GTOWaveFunc[i]
        ne *= ne
        ns += ne
     
    return ns

