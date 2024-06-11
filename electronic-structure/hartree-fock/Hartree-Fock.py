from abc import abstractmethod, ABC


class HF(ABC):
    def __init__(self,Z,basis):
        self.Z = Z
        self.basis = basis

    def findMatricies(self):
        pass

class RHF(HF):
    def __init__(self,Z,basis):
        HF.__init__(self,Z,basis)

class UHF(HF):
    def __init__(self,Z,basis):
        HF.__init__(self,Z,basis)