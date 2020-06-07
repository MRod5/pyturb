"""
Gas:
----

Class for a generic gas. Object Gas has Rg and Ru properties and methods to calculated
specific heat capacities (cp, cv) and heat capacity ratio.

MRodriguez. 2020

"""

from abc import abstractmethod


class Gas(object):
    """
    Class for a generic gas.

    + Properties:
        gas_species: Species selected. May be a pure substance or any of the molecules and mixes
                     considered in "NASA Glenn Coefficients for Calculating Thermodynamic
                     Properties of Individual Species".
        Ru: Ideal gas law universal constant  [J/mol/K]
        Rg: Individual gas constant [J/kg/K]
    + Methods:
        cp: Heat capacity at constant pressure [J/kg/K]
        cv: Heat capacity at constant volume [J/kg/K]
        gamma: Heat capacity ratio [-]

    """

    @property
    @abstractmethod
    def gas_species(self):
        """
        Gets the Name of the gas species selected. May be a pure substance or any of the 
        molecules and mixes considered in "NASA Glenn Coefficients for Calculating Thermodynamic
        Properties of Individual Species".
        """
        raise NotImplementedError

        
    @property
    @abstractmethod
    def Ru(self):
        """
        Get the Ideal Gas Law constant Ru [J/mol/K]
        """
        return NotImplementedError


    @property
    @abstractmethod
    def Rg(self):
        """
        Get the Individual Gas constant Rg =  Ru/Mg [J/kg/K]
        """
        raise NotImplementedError


    @abstractmethod
    def cp(self, temperature):
        """
        Heat capacity ratio at constant pressure [J/kg/K]
        """
        raise NotImplementedError


    @abstractmethod
    def cv(self, temperature):
        """
        Heat capacity ratio at constant volume [J/kg/K]
        """
        raise NotImplementedError


    @abstractmethod
    def gamma(self, temperature):
        """
        Heat capacity ratio cp/cv [-]
        """
        raise NotImplementedError

