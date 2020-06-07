"""
Semiperfect Ideal Gas:
----------------------

Semiperfect gas approach with Ideal Gas: cp, cv are calculated as a function of temperature,
and the Ideal Gas Law (pv=RgT) and the Mayer Equation (cp-cv=Rg) are applicable:
    u = cv(T)·T
    h = cp(T)·T
    pv = Rg·T
    h = u + pv --> cp(T)·T = cv(T)·T + Rg·T --> Rg = cp(T) - cv(T)

MRodriguez. 2020

"""

import numpy as np
import pyturb.utils.constants as cts
from pyturb.gas_models.gas import Gas
from pyturb.gas_models.thermo_properties import ThermoProperties


class SemiperfectIdealGas(Gas):
    """
    SemiperfectIdealGas class:
    ----------------------

    Implements a Gas object:
        +gas_species: Species selected. May be a pure substance or any of the molecules and mixes
                      considered in "NASA Glenn Coefficients for Calculating Thermodynamic
                      Properties of Individual Species".
        + Ru: Ideal gas law universal constant  [J/mol/K]
        + Rg: Individual gas constant [J/kg/K]
        + cp: Heat capacity at constant pressure [J/kg/K]
        + cv: Heat capacity at constant volume [J/kg/K]
        + gamma: Heat capacity ratio [-]
    
    When initialized, ThermoProperties is called with the selected gas species to retrieve
    the properties of the gas: Mg, deltaHf_ef, deltaHf_0K and coefficients for calculating cp(T).

    Heat capacity at constant pressue is obtained as a polynomial with 7 coefficients:
        cp/Rg = a1*T**(-2) + a2*T**(-1) + a3 + a4*T**(1) + a5*T**(2) + a6*T**(3) + a7*T**(4)

    Heat capacity at constant volume and heat capacity ratio are obtained as an Ideal Gas:
    cv = cp - Rg, gamma=cp/cv.
    
    """

    def __init__(self, species=''):
        """
        Gets the thermodynamic properties and the gas species name.
        """

        # Check if the gas species is valid
        if species=='':
            raise ValueError("Gas species not defined.")
        if type(species) is not str:
            raise TypeError("Gas species must be a string")


        # Get thermodynamic properties
        self.thermo_prop = ThermoProperties(species)

        if self.thermo_prop.Mg == None:
            # If gas species is not valid:
            raise ValueError("Gas unavailable: {}".format(species))
        else:
            self.__gas_species = species
        return
        
    
    @property
    def gas_species(self):
        """
        Gets the Name of the gas species selected. May be a pure substance or any of the 
        molecules and mixes considered in "NASA Glenn Coefficients for Calculating Thermodynamic
        Properties of Individual Species".
        """
        return self.__gas_species
        

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
        Rg = self.Ru/self.thermo_prop.Mg*1e3
        return Rg
        
    
    def cp(self, temperature):
        """
        Heat capacity ratio at constant pressure [J/kg/K].

        As a semiperfect gas, cp is a function of the temperature. It is calculated as a
        7 coefficients polynomial:
            cp/Rg = a1*T**(-2) + a2*T**(-1) + a3 + a4*T**(1) + a5*T**(2) + a6*T**(3) + a7*T**(4)
        """

        # Check if temperature is out of range for the current gas species:
        tmin = np.min(self.thermo_prop.temp_range)
        tmax = np.max(self.thermo_prop.temp_range)

        if tmin>temperature or tmax<temperature:
            raise ValueError("Gas temperature ({0}K) out of implemented limits [{1},{2}]K".format(temperature, tmin, tmax))
        else:
            for (temp_range, coeffs) in zip(self.thermo_prop.temp_range, self.thermo_prop.coefficients):
                # For all the temperature intervals of the gas species, pick the correct one:
                # min_T_interval <= temperature max_T_interval

                if temp_range[0]<=temperature<temp_range[1]:
                    # Calculate 7 terms polynomial (non-dimensional):
                    temp_poly = np.array([temperature**(-2), temperature**(-1), 1, temperature, temperature**(2), temperature**(3), temperature**(4)])
                    cp_ = np.dot(coeffs, temp_poly)
                    
                    # Dimensional cp_ at current temperature:
                    cp_ = cp_ * self.Rg
        return cp_
        
    
    def cv(self, temperature):
        """
        Heat capacity ratio at constant volume [J/kg/K]

        As a semiperfect gas, cv is a function of tempeature. cv is calculated with the 
        Mayer equation: cv(T) = cp(T) - Rg (ideal gas).
        """
        cv_ = self.cp(temperature) - self.Rg
        return cv_
    
        
    def gamma(self, temperature):
        """
        Heat capacity ratio cp/cv [-].

        As a perfect gas, gamma is considered invariant with temperature. gamma is calculated
        as gamma = cp/cv.
        """
        
        gamma_ = self.cp(temperature)/self.cv(temperature)
        return gamma_
        
        