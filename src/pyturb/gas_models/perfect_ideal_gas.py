"""
Perfect Ideal Gas:
------------------

Perfect gas approach with Ideal Gas: cp, cv are constant and the Ideal Gas Law (pv=RgT) and
the Mayer Equation (cp-cv=Rg) are applicable:

    u = cv·T
    h = cp·T
    pv = Rg·T
    h = u + pv --> cp·T = cv·T + Rg·T --> Rg = cp - cv

MRodriguez. 2020

"""

import numpy as np
import pyturb.utils.constants as cts
from pyturb.gas_models.gas import Gas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas


class PerfectIdealGas(Gas):
    """
    PerfectIdealGas class:
    ----------------------

    Implements a Gas object:
        +gas_species: Species selected. May be a pure substance or any of the molecules and mixes
                      considered in "NASA Glenn Coefficients for Calculating Thermodynamic
                      Properties of Individual Species".
        + Ru: Ideal gas law universal constant  [J/mol/K]
        + Rg: Individual gas constant [J/kg/K]
        + cp: Constant heat capacity at constant pressure [J/kg/K]
        + cv: Constant heat capacity at constant volume [J/kg/K]
        + gamma: Heat capacity ratio [-]
        + cp_molar: Constant molar heat capacity at constant pressure [J/kg/K]
        + cv_molar: Constant molar heat capacity at constant volume [J/kg/K]
        + h0_molar: Assigned molar enthalpy [J/mol] as h0(T) = deltaHf(T_ref) + (H0(T) - h0(T_ref))
        + h0: Assigned enthalpy [J/kg] as h0(T) = deltaHf(T_ref) + (H0(T) - h0(T_ref))
    
    When initialized, SemiperfectIdealGas class is called to retrieve the gas species.
    The heat capacity at constant pressure is calculated as SemiperfectIdealGas.cp(T_ref)
    where T_ref is the reference temperatue for chemical processes.

    Heat capacity at constant volume and heat capacity ratio are obtained as an Ideal Gas:
    cv = cp - Rg, gamma=cp/cv.
    
    """

    def __init__(self, species=''):
        """
        Initializes a SemiperfectIdealGas object with the gas species selected and stores the
        gas species and the reference temperature.
        """

        # Check if the gas species is valid
        if species=='':
            raise ValueError("Gas species not defined.")
        if type(species) is not str:
            raise TypeError("Gas species must be a string")

        # SemiperfectIdealGas object:
        self.__from_semiperfect_gas = SemiperfectIdealGas(species)
        
        # Store the thermodynamic properties of the gas
        self.thermo_prop = self.__from_semiperfect_gas.thermo_prop

        # Reference temperature and gas species
        self.T_ref = cts.T_ref
        self._gas_species = species
        
        return
        
    
    @property
    def gas_species(self):
        """
        Gets the Name of the gas species selected. May be a pure substance or any of the 
        molecules and mixes considered in "NASA Glenn Coefficients for Calculating Thermodynamic
        Properties of Individual Species".
        """
        
        return self._gas_species
        

    @property
    def Ru(self):
        """
        Get the Ideal Gas Law constant Ru [J/mol/K]
        """
        
        Ru = cts.Ru
        return Ru
        
        
    @property
    def Rg(self):
        """
        Get the Individual Gas constant Rg =  Ru/Mg [J/kg/K]
        """
        
        Rg = self.Ru/self.Mg*1e3
        return Rg


    @property
    def Mg(self):
        """
        Get the molecular mass Mg [g/mol]
        """
        
        return self.thermo_prop.Mg
        
        
    def cp_dimensionless(self, temperature=None):
        """
        Dimensionless heat capacity ratio at constant pressure [-].

        As a perfect gas, cp is considered invariant with temperature. It is calculated as a
        semi-perfect gas (which means cp(T)) at T_ref temperature for any temperature.
        """
        
        cp_ = self.__from_semiperfect_gas.cp_dimensionless(self.T_ref)
        return cp_


    def gamma(self, temperature=None):
        """
        Heat capacity ratio cp/cv [-].

        As an ideal gas, gamma is considered invariant with temperature. gamma is calculated
        as gamma = cp/cv.
        """
        
        gamma_ = self.cp(self.T_ref)/self.cv(self.T_ref)
        return gamma_


    def cp(self, temperature=None):
        """
        Heat capacity ratio at constant pressure [J/kg/K].

        As a perfect gas, cp is considered invariant with temperature. It is calculated as a
        semi-perfect gas (which means cp(T)) at T_ref temperature for any temperature.
        """
        
        cp_ = self.cp_dimensionless()*self.Rg
        return cp_
     

    def cv(self, temperature=None):
        """
        Heat capacity ratio at constant volume [J/kg/K]

        As an ideal gas, cv is invariant with tempeature. cv is calculated with the 
        Mayer equation: cv = cp - Rg.
        """
        
        cv_ = self.cp() - self.Rg
        return cv_
    
        
    def cp_molar(self, temperature=None):
        """
        Molar heat capacity ratio at constant pressure [J/mol/K].

        As a perfect gas, cp is considered invariant with temperature. It is calculated as a
        semi-perfect gas (which means cp(T)) at T_ref temperature for any temperature.
        """
        
        cp_ = self.cp_dimensionless()*self.Ru
        return cp_
     

    def cv_molar(self, temperature=None):
        """
        Molar heat capacity ratio at constant volume [J/mol/K]

        As an ideal gas, cv is invariant with tempeature. cv is calculated with the 
        Mayer equation: cv = cp - Ru.
        """
        
        cv_ = self.cp_molar() - self.Ru
        return cv_


    def h0_dimensionless(self, temperature=None):
        """
        Dimensionless assigned enthalpy:
            h0(T)/R = deltaHf(T_ref)/R + (H0(T) - h0(T_ref))/R [-].
        Dimensionless, semiperfect gas specific heat capacity at constant pressure 
        is used.
        """
        
        h0_ = self.cp_dimensionless()

        return h0_
        
        
    def h0(self, temperature):
        """
        h0(T) = deltaHf(T_ref) + (H0(T) - h0(T_ref)) = int(cp(T)*dT). [J/kg]
        Perfect gas specific heat capacity at constant pressure is used.
        """
        
        h_ = self.h0_dimensionless(temperature) * self.Rg * temperature
                    
        return h_
    
    
    def h0_molar(self, temperature):
        """
        h0(T) = deltaHf(T_ref) + (H0(T) - h0(T_ref)) = int(cp(T)*dT). [J/mol]
        Perfect gas specific heat capacity at constant pressure is used.
        """
        
        h_ = self.h0_dimensionless(temperature) * self.Ru * temperature
                    
        return h_