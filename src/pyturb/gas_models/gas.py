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
        Ng: Number of moles of the gas [mol]
        mg: Mass quantity of the gas [kg]
        Mg: Molecular mass of the gas [g/mol]
        Ru: Ideal gas law universal constant  [J/mol/K]
        Rg: Individual gas constant [J/kg/K]
    + Methods:
        cp: Heat capacity at constant pressure [J/kg/K]
        cv: Heat capacity at constant volume [J/kg/K]
        cp_molar: Molar heat capacity at constant pressure [J/kg/K]
        cv_molar: Molar heat capacity at constant volume [J/kg/K]
        gamma: Heat capacity ratio [-]
        h0: Assigned enthalpy [J/kg]
        h0_molar: Assigned molar enthalpy [J/mol]

    """


    @property
    @abstractmethod
    def gas_species(self):
        """
        Gets the Name of the gas species selected. May be a pure substance or any of the 
        molecules and mixes considered in "NASA Glenn Coefficients for Calculating Thermodynamic
        Properties of Individual Species".
        """
        return NotImplementedError


    @property
    @abstractmethod
    def Ng(self):
        """
        Get the number of moles of gas [mol]
        """
        return NotImplementedError
    

    @property
    @abstractmethod
    def mg(self):
        """
        Get the mass quantity of gas [kg]
        """
        return NotImplementedError


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


    @property
    @abstractmethod
    def Mg(self):
        """
        Molecular mass of the gas. [g/mol]
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


    @abstractmethod
    def cp_molar(self, temperature):
        """
        Molar heat capacity ratio at constant pressure [J/mol/K]
        """
        raise NotImplementedError


    @abstractmethod
    def cv_molar(self, temperature):
        """
        Molar heat capacity ratio at constant volume [J/mol/K]
        """
        raise NotImplementedError
        
    
    @abstractmethod
    def h0(self, temperature):
        """
        Assigned enthalpy h0(T) = deltaHf(T_ref) + (H0(T) - h0(T_ref)) [J/kg]
        """
        raise NotImplementedError
    
    
    @abstractmethod
    def h0_molar(self, temperature):
        """
        Assigned molar enthalpy h0(T) = deltaHf(T_ref) + (H0(T) - h0(T_ref)) [J/mol]
        """
        raise NotImplementedError