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
        + cp_molar: Molar heat capacity at constant pressure [J/mol/K]
        + cv_molar: Molar heat capacity at constant volume [J/mol/K]
        + h0_molar: Assigned molar enthalpy [J/mol] as h0(T) = deltaHf(T_ref) + (H0(T) - h0(T_ref))
        + h0: Assigned enthalpy [J/kg] as h0(T) = deltaHf(T_ref) + (H0(T) - h0(T_ref))
    
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
        

    @property
    def Mg(self):
        """
        Get the molecular mass Mg [g/mol]
        """
        
        return self.thermo_prop.Mg

    
    def cp_dimensionless(self, temperature):
        """
        Dimensionless heat capacity ratio at constant pressure [-].

        As a semiperfect gas, cp is a function of the temperature. It is calculated as a
        7 coefficients polynomial:
            cp/Rg = a1*T**(-2) + a2*T**(-1) + a3 + a4*T**(1) + a5*T**(2) + a6*T**(3) + a7*T**(4)
        """

        # Check if temperature is out of range for the current gas species:
        tmin = np.min(self.thermo_prop.temp_range)
        tmax = np.max(self.thermo_prop.temp_range)

        if np.any(tmin>temperature) or np.any(tmax<temperature):
            raise ValueError("Gas temperature ({0}K) out of implemented limits [{1},{2}]K".format(temperature, tmin, tmax))
        else:
            cp_ = np.zeros_like(temperature)
            for (temp_range, coeffs) in zip(self.thermo_prop.temp_range, self.thermo_prop.coefficients):
                # For all the temperature intervals of the gas species, pick the correct one:
                # min_T_interval <= temperature max_T_interval

                mask_temp = np.logical_and(temp_range[0]<=temperature, temp_range[1]>= temperature)
                temp_poly = np.array([temperature**(-2), temperature**(-1), 1, temperature, temperature**(2), temperature**(3), temperature**(4)], dtype="object")
                cp_ = cp_ + np.dot(coeffs, temp_poly) * mask_temp
            
        return cp_


    def gamma(self, temperature):
        """
        Heat capacity ratio cp/cv [-].

        As a perfect gas, gamma is considered invariant with temperature. gamma is calculated
        as gamma = cp/cv.
        """
        
        gamma_ = self.cp(temperature)/self.cv(temperature)
        return gamma_


    def cp(self, temperature):
        """
        Heat capacity ratio at constant pressure [J/kg/K].

        As a semiperfect gas, cp is a function of the temperature. It is calculated as a
        7 coefficients polynomial:
            cp/Rg = a1*T**(-2) + a2*T**(-1) + a3 + a4*T**(1) + a5*T**(2) + a6*T**(3) + a7*T**(4)
        """

        cp_ = self.cp_dimensionless(temperature) * self.Rg
        return cp_


    def cp_molar(self, temperature):
        """
        Molar heat capacity ratio at constant pressure [J/mol/K].

        As a semiperfect gas, cp is a function of the temperature. It is calculated as a
        7 coefficients polynomial:
            cp/Rg = a1*T**(-2) + a2*T**(-1) + a3 + a4*T**(1) + a5*T**(2) + a6*T**(3) + a7*T**(4)
        """

        cp_ = self.cp_dimensionless(temperature) * self.Ru
        return cp_

    
    def cv(self, temperature):
        """
        Heat capacity ratio at constant volume [J/kg/K]

        As a semiperfect gas, cv is a function of tempeature. cv is calculated with the 
        Mayer equation: cv(T) = cp(T) - Ru (ideal gas).
        """
        cv_ = self.cp(temperature) - self.Rg
        return cv_


    def cv_molar(self, temperature):
        """
        Molar heat capacity ratio at constant volume [J/kg/K]

        As a semiperfect gas, cv is a function of tempeature. cv is calculated with the 
        Mayer equation: cv(T) = cp(T) - Rg (ideal gas).
        """
        cv_ = self.cp_molar(temperature) - self.Ru
        return cv_


    def h0_dimensionless(self, temperature):
        """
        Dimensionless assigned enthalpy:
            h0(T)/R = deltaHf(T_ref)/R + (H0(T) - h0(T_ref))/R [-].
        Dimensionless, semiperfect gas specific heat capacity at constant pressure 
        is used.
        """
        
        # Check if temperature is out of range for the current gas species:
        tmin = np.min(self.thermo_prop.temp_range)
        tmax = np.max(self.thermo_prop.temp_range)

        if tmin>temperature or tmax<temperature:
            raise ValueError("Gas temperature ({0}K) out of implemented limits [{1},{2}]K".format(temperature, tmin, tmax))
        else:
            for (temp_range, coeffs, integ_cts) in zip(self.thermo_prop.temp_range, self.thermo_prop.coefficients, self.thermo_prop.integration_cts):
                # For all the temperature intervals of the gas species, pick the correct one:
                # min_T_interval <= temperature max_T_interval

                if temp_range[0]<=temperature<temp_range[1]:
                    # Calculate 8 terms polynomial (non-dimensional):
                    temp_poly = np.array([-temperature**(-2), np.log(temperature)/temperature, 1, temperature**(1)/2, temperature**(2)/3, temperature**(3)/4, temperature**(4)/5])
                    h0_ = np.dot(coeffs, temp_poly) + integ_cts[0] / temperature

            return h0_
        
        
    def h0(self, temperature):
        """
        Assigned enthalpy:
            h0(T) = deltaHf(T_ref) + (H0(T) - h0(T_ref)) = int(cp(T)*dT). [J/kg]
       
        """
            
        h0_ = self.h0_dimensionless(temperature) * self.Rg * temperature
            
        return h0_
            
        
    def h0_molar(self, temperature):
        """
        Assigned molar enthalpy:
            h0(T) = deltaHf(T_ref) + (H0(T) - h0(T_ref)) = int(cp(T)*dT). [J/mol]

        """
            
        h0_ = self.h0_dimensionless(temperature) * self.Ru * temperature
            
        return h0_
    
    