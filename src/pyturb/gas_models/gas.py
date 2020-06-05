"""
"""

from abc import abstractmethod


class Gas(object):
    """
    """

    @property
    @abstractmethod
    def gas_species(self):
        """
        """
        raise NotImplementedError

        
    @property
    @abstractmethod
    def Ru(self):
        """
        Get the Ideal Gas law constant Rg =  Ru/Mg
        """
        return NotImplementedError


    @property
    @abstractmethod
    def Rg(self):
        """
        """
        raise NotImplementedError


    @abstractmethod
    def cp(self, temperature):
        """
        """
        raise NotImplementedError


    @abstractmethod
    def cv(self, temperature):
        """
        """
        raise NotImplementedError


    @abstractmethod
    def gamma(self, temperature):
        """
        """
        raise NotImplementedError

