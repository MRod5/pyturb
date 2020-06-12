"""
"""

from abc import abstractmethod
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas



class ControlVolume(object):
    """
    """
    def __init__(self, fluid):
        if not(isinstance(fluid, PerfectIdealGas) or isinstance(fluid, SemiperfectIdealGas)):
            raise TypeError("PerfectIdealGas, SemiperfectIdealGas")
        self.fluid = fluid
        return

    @abstractmethod
    def het(self):
        """
        """
        raise NotImplementedError


    @abstractmethod
    def hst(self):
        """
        """
        raise NotImplementedError

    
    @property
    @abstractmethod
    def Tet(self):
        """
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def Tst(self):
        """
        """
        raise NotImplementedError
