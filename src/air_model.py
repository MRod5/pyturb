#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
air_model:
---------

This file contains two approaches for the air: one assuming it is a thermally perfect gas, 
and one with polynomials for specific heat and heat capacity ratio (both at constant pressure and volume).

Content:
--------
    + Class Air:
        + R_univ : float. Universal Gas Constant [J/mol/K]
        + M_O: float. Molecular mass of the oxygen [g/mol] (monoatomic)
        + M_N_ float. Molecular mass of the nitrogen [g/mol] (monoatomic)
        + M_air: float. Molecular mass of the air [g/mol]
        + R_air: float. Specific air constant [J/kg/K]
        + cp_air: Specific heat capacity to be implemented ar RealGas or IdealGas
        + cv_air: Specific heat capacity to be implemented ar RealGas or IdealGas
        + gamma_air: Specific heat capacity to be implemented ar RealGas or IdealGas
    
    + class IdealGas:
        + cp [J/kg/K]
        + cv [J/kg/K]
        + gamma_air [non-dimensional]
    
    + class RealGas:
        + cp_lp: function. Low pressure cp for air
        + cp_naca1135: empirical formula for cp as a function of temperature [J/kg/K]
        + cp_nasa4513: Nasa technical memorandum, thermodynamic data coefficients [J/kg/K]

Dependencies:
-------------
    + numpy
    + kelvin2rankine, rankine2kelvin and ft2m from units
    + interp1d from scipy.interpolate



Part of pyturb. Tests pending.

April 2017

MRodriguez
"""

import numpy as np
from abc import abstractmethod
from units import kelvin2rankine, rankine2kelvin, ft2m
from scipy.interpolate import interp1d


class Air(object):
    """
    Air class:
    ----------
    
    Molar masses, air composition and abstract methods for cp_air, cv_air and gamma_air.
    
    Contents:
    ---------
        + R_univ : float. Universal Gas Constant [J/mol/K]
        + M_O: float. Molecular mass of the oxygen [g/mol] (monoatomic)
        + M_N_ float. Molecular mass of the nitrogen [g/mol] (monoatomic)
        + M_air: float. Molecular mass of the air [g/mol]
        + R_air: float. Specific air constant [J/kg/K]
        + cp_air: Specific heat capacity to be implemented ar RealGas or IdealGas
        + cv_air: Specific heat capacity to be implemented ar RealGas or IdealGas
        + gamma_air: Specific heat capacity to be implemented ar RealGas or IdealGas
    
    """

    def __init__(self):
        """
        Initializes Air class. Universal gas constant, molar masses for oxygen, nitrogen and air and specific air 
        constant are defined here.
        """
        # Universal gas constant:
        self.R_univ = 8.3144598  # J/mol/K

        # Molecular masses:
        # Define molecular masses for oxygen and nitrogen and calculate molecular mass of the air
        # Air is assumed to be composed of 22% oxygen and 78% nitrogen:

        self.M_O = 15.9994  # Da
        self.M_N = 14.0067  # Da
        # self.M_Ar = 39.94   # Da

        # self.M_air = 0.2095 * 2 * self.M_O + 0.78 * 2 * self.M_N + 0.0093 * self.M_Ar # Da

        self.M_air = 0.22*2*self.M_O + 0.78*2*self.M_N

        # Air specific constant:
        self.R_air = self.R_univ / (self.M_air * 1e-3)  # IS Air specific constant [J/kg/K]

        return

    @abstractmethod
    def cp_air(self, temp):
        """
        Specific heat capacity at constant pressure
        """
        pass

    @abstractmethod
    def cv_air(self, temp):
        """
        Specific heat capacity at constant volume
        """
        pass

    @abstractmethod
    def gamma_air(self, temp):
        """
        Heat capacity ratio
        """
        pass


class IdealGas(Air):
    """
    Heat capacity ratio and Specific heat capacity at constant pressure and constant volume for air assuming
    ideal gas.
    
    Inherits Air superclass. 
    
    Contents:
    ---------
        cp_air: float. Specific heat capacity at constant pressure for air as an ideal gas [J/kg/K]
        cv_air: float. Specific heat capacity at constant volume for air as an ideal gas [J/kg/K]
        gamma_air: Specific heat ratio for air as an ideal gas [non-dimensional]
    """

    def __init__(self):
        """
        Initializes IdealGas properties.
        """
        super().__init__()
        self.cp_air_ = 1004  # J/kg/K
        self.cv_air_ = self.cp_air_ - self.R_air  # Mayer equation with the air constant [J/kg/K]
        self.gamma_air_ = self.cp_air_ / self.cv_air_  # non-dimensional

    def cp_air(self, temp=None):
        """
        Specific heat capacity at constant pressure
        
        Inputs:
        -------
            + temp: Current temperature
        Outputs:
        --------
            + cp_air: float. Constant cp for air (perfect gas)
        """
        return self.cp_air_

    def cv_air(self, temp=None):
        """
        Specific heat capacity at constant volume

        Inputs:
        -------
            + temp: Current temperature
        Outputs:
        --------
            + cv_air: float. Constant cv for air (perfect gas)
        """
        return self.cv_air_

    def gamma_air(self, temp=None):
        """
        Heat capacity ratio

        Inputs:
        -------
            + temp: Current temperature
        Outputs:
        --------
            + gamma: float. Constant gamma for air (perfect gas)
        """
        return self.gamma_air_


class RealGas(Air):
    """
    Heat capacity ratio and Specific heat capacity at constant pressure and constant volume for air as a function of
    temperature.
    
    RealGas(cp_option, gamma_option):
    ---------------------------------
        + cp_option: 'low-pressure', 'nasa' or 'naca'. Cp function selector. Default is 'nasa'
        + gamma_option: 'nasa' or 'standard'. Gamma_air function selector. Default is 'standard'
    
    
    Contents:
    ---------
        + cp_air: calculates the cp of air with the function selected at cp_option   
        + cv_air: calculates de cv of air with Mayer equation (cp_air calculated with cp_option)
        + gamma_air: heat capacity ratio obtained with gamma_option function
        + cp_lp: low pressure specific heat at constant pressure as a 3rd degree polynomial 
        + cp_naca1135: NACA empirical equation for specific heat at constant pressure
        + cp_nasa4513: NASA 4th degree polynomial for specific heat at constant pressure
        + gamma_air_naca: NACA piece-wise function for heat capacity ratio, linearly interpolated depending 
                on input temperature
        + gamma_air_standard: Piece-wise heat capacity ratio as a function of temperature
    """

    def __init__(self, cp_option='nasa', gamma_option='standard'):
        super().__init__()

        # Dictionary for cp functions:
        self.cp_air_functions = {'low-pressure': self.cp_lp, 'nasa': self.cp_nasa4513, 'naca': self.cp_naca1135}

        # Dictionary for heat capacity ratio functions:
        self.gamma_air_functions = {'naca': self.heat_capacity_ratio_naca, 'standard': self.heat_capacity_ratio_std}

        # Check selected options:
        # gamma_air:
        if type(gamma_option) == str:
            if gamma_option.lower() in self.gamma_air_functions.keys():
                # Selected option for calculating heat capacity ratio:
                self.gamma_option = gamma_option
            else:
                print('Wrong gamma_air function {}. Selecting standard gamma_air instead'.format(gamma_option))
                self.gamma_option = 'standard'
        else:
            print('Wrong gamma_air function {}. Selecting standard gamma_air instead'.format(gamma_option))
            self.gamma_option = 'standard'

        # cp_air:
        if type(cp_option) == str:
            if cp_option.lower() in self.cp_air_functions.keys():
                # Selected option for calculating cp, cv:
                self.cp_option = cp_option
            else:
                print('Wrong cp_air function {}. Selecting NASA4513 instead'.format(cp_option))
                self.cp_option = 'nasa'
        else:
            print('Wrong cp_air function {}. Selecting NASA4513 instead'.format(cp_option))
            self.cp_option = 'nasa'

        # Class parameters:
        self.theta = rankine2kelvin(5500)  # Reference temperature for empirical cp function
        self.cp_air_perfect = 6006  # Perfect cp_air for empirical cp function [ft**2/s**2/R]
        self.cp_air_perfect *= ft2m ** 2 * kelvin2rankine(1)   # Converted to SI: [J/kg/K]

        self.gamma_air_perfect = 7 / 5  # Perfect specific heat ratio
        self.gamma_air_comb = 1.3  # High temperature characteristic specific heat ratio
        self.gamma_air_comp = 1.38  # Compression zone characteristic specific heat ratio
        self.gamma_air_exp = 1.33  # Expansion zone characteristic specific heat ratio

    def cp_air(self, temp):
        """
        Specific heat capacity at constant pressure selector:
            + cp_lp
            + cp_
        
        Inputs:
        -------
            + temp: Temperature at which the cp of air will be calculated [K]
        
        Outputs:
        --------
            + cp_air: Specific heat at constant pressure [J/kg/K]
        """
        return self.cp_air_functions[self.cp_option](temp)

    def cv_air(self, temp):
        """
        Specific heat capacity at constant volume obtained with the Mayer Equation (Cp - R). cp is obtained with cp_air.

        Inputs:
        -------
            + temp: Temperature at which the cp of air will be calculated [K]

        Outputs:
        --------
            + cv_air: Specific heat at constant volume [J/kg/K]
        """
        return self.cp_air(temp) - self.R_air

    def gamma_air(self, temp):
        """
        Heat capacity ratio selector. 

        Inputs:
        -------
            + temp: Temperature at which the cp of air will be calculated [K]

        Outputs:
        --------
            + gamma_air: Heat capacity ratio [non-dimensional]
        """
        return self.gamma_air_functions[self.gamma_option](temp)

    def cp_lp(self, temp):
        """
        Specific heat capacity at constant pressure for air at low pressure (vacuum pressures) as a function of the 
        temperature by means of a 3rd degree polynomial.
        
        + Inputs:
        ---------
            temp: float. Temperature to calculate the cp [K]
            
        + Outputs:
        ----------
            cp_air: float. Heat capacity at constant pressure [J/kg/K]
 
        """
        # Note tat the coefficients are expressed in mol. Conversion to kg is made before the return with the molecular
        # of the air calculated above.

        # Coefficients:
        cp0 = 28.11  # J/mol/K
        cp1 = 0.1967e-2  # J/mol/K
        cp2 = 0.4802e-4  # J/mol/K
        cp3 = -1.966e-9  # J/mol/K

        # cp:
        cp_air = (cp0 + cp1 * temp + cp2 * temp ** 2 + cp3 * temp ** 3) / (self.M_air * 1e-3)  # J/mol/K
        return cp_air

    def cp_naca1135(self, temp):
        """
        
        Specific heat capacity at constant pressure for air by means of an empirical equation. See NACA 1135 report.
        
        The empirical equation is valid from 400ºR to 5000ºR.
        
        + Inputs:
        ---------
            temp: float. Temperature to calculate the cp [K]
            
        + Outputs:
        ----------
            cp_air: float. Empirical air capacity at constant pressure [J/kg/K]
        
        """

        if not rankine2kelvin(385) <= temp <= rankine2kelvin(5000):
            print("Temperature range not yet implemented")
            return

        # Temperature coefficient for empirical equation
        temp_rel = self.theta / temp

        # Empirical cp:
        cp_air = self.cp_air_perfect * (1 + temp_rel ** 2 * np.exp(temp_rel) /
                                        (np.exp(temp_rel)) ** 2 * (self.gamma_air_perfect - 1) / self.gamma_air_perfect)

        return cp_air

    def cp_nasa4513(self, temp):
        """
        
        Specific heat capacity at constant pressure for air as a function of the local temperature by means of a 4th
        degree polynomial. See NASA Technical Memorandum 4513 for further details.
        
        Two sets of coefficients are applied depending on the input temperature. From 200 to 1000K and from 1000K to
        6000K respectively. The polynomial is non-dimensional, dimension are obtained with the specific constant of
        the air.
        
        + Inputs:
        ---------
            temp: float. Temperature to calculate the cp [K]
            
        + Outputs:
        ----------
            cp_air: float. Air capacity at constant pressure as a function of the temperature [J/kg/K]
        """

        # Apply coefficients depending of the temperature:
        if 1000 <= temp <= 6000:
            a0_air = 3.08792717E+00  # Independent coefficient
            a1_air = 1.24597184E-03  # 1st order [1/K]
            a2_air = -4.23718945E-07  # 2nd order [1/K**2]
            a3_air = 6.74774789E-11  # 3rd order [1/K**3]
            a4_air = -3.97076972E-15  # 4th order [1/K**4]
        elif 200 <= temp < 1000:
            a0_air = 3.56839620E+00  # Independent coefficient
            a1_air = -6.78729429E-04  # 1st order [1/K**4]
            a2_air = 1.55371476E-06  # 2nd order [1/K**4]
            a3_air = -3.29937060E-12  # 3rd order [1/K**4]
            a4_air = -4.66395387E-13  # 4th order [1/K**4]
        else:
            print("Temperature range not yet implemented")
            return

        # Apply polynomial:
        cp_air = a0_air + a1_air * temp + a2_air * temp ** 2 + a3_air * temp ** 3 + a4_air * temp ** 4

        # Get dimensions (IS)
        return cp_air * self.R_air  # [J/kg/K]

        # N2 REF ELEMENT    G 8/02N  2.   0.   0.   0.G   200.000  6000.000 1000.        1
        #  2.95257637E+00 1.39690040E-03-4.92631603E-07 7.86010195E-11-4.60755204E-15    2
        # -9.23948688E+02 5.87188762E+00 3.53100528E+00-1.23660988E-04-5.02999433E-07    3
        #  2.43530612E-09-1.40881235E-12-1.04697628E+03 2.96747038E+00 0.00000000E+00    4
        #
        # O2 REF ELEMENT    RUS 89O   2    0    0    0G   200.000  6000.000 1000.        1
        #  3.66096083E+00 6.56365523E-04-1.41149485E-07 2.05797658E-11-1.29913248E-15    2
        # -1.21597725E+03 3.41536184E+00 3.78245636E+00-2.99673415E-03 9.84730200E-06    3
        # -9.68129508E-09 3.24372836E-12-1.06394356E+03 3.65767573E+00 0.00000000E+00    4
        #
        # AIR               L 9/95  WARNING!         0G   200.000  6000.000 1000.        1
        #  3.08792717E+00 1.24597184E-03-4.23718945E-07 6.74774789E-11-3.97076972E-15    2
        # -9.95262755E+02 5.95960930E+00 3.56839620E+00-6.78729429E-04 1.55371476E-06    3
        # -3.29937060E-12-4.66395387E-13-1.06234659E+03 3.71582965E+00-1.50965000E+01    4

    @staticmethod
    def heat_capacity_ratio_naca(temp):
        """
        Static method.
        Heat capacity ratio as a function of the temperature. See NACA1135 report for further details.
          
        The heat capacity ratio in this function is valid from 400ºR to 5000ºR. At that range, any input temperature 
        will be used to interpolate the value of gamma_air.
        
        + Inputs:
        ---------
            temp: float. Temperature to calculate the heat capacity ratio [K]
            
        + Outputs:
        ----------
            gamma_air: array. Empirical air heat capacity ratio 
            
        + Dependencies:
        ---------------
            interp1d from scipy.interpolate
        """

        if not rankine2kelvin(380) <= temp <= rankine2kelvin(5000):
            print("Temperature range not yet implemented")
            return

        # Known heat capacity ratios
        gamma_air_fixed = np.array([1.4, 1.399, 1.396, 1.392, 1.387, 1.381, 1.375, 1.368,
                                    1.361, 1.355, 1.349, 1.344, 1.339, 1.335, 1.331, 1.328,
                                    1.322, 1.317, 1.313, 1.309, 1.306, 1.301, 1.298, 1.296, 1.294])

        # Temperatures of the known heat capacity ratios
        temp_fixed = np.concatenate([np.linspace(500, 2000, 16), np.arange(2200, 3000, 200),
                                     np.arange(3000, 5500, 500)])

        # Linear interpolation. If temperature is outside bounds, the interpolation wont work (shouldn't happen,
        # temperature check already passed)
        gamma_air_function = interp1d(temp_fixed, gamma_air_fixed, kind='linear')

        # If temperature is below 277.7K, gamma_air is fixed to 1.4:
        if kelvin2rankine(temp) < temp_fixed[0]:
            return gamma_air_fixed[0]
        gamma_air = gamma_air_function(kelvin2rankine(temp))

        return float(gamma_air)

    def heat_capacity_ratio_std(self, temp):
        """
        Piece-wise heat capacity ratio as a function of the temperature. 
          
        The heat capacity ratio is made constant during the temperature interval where it is defined:
        
        Until 600K: gamma_air_perfect, from 600 to 900K gamma_air_comp, from 900K to 1400K gamma_air_exp and 
        from 1400K to 2000K gamma_air_comb.
        
        + Inputs:
        ---------
            temp: float. Temperature to calculate the heat capacity ratio [K]
            
        + Outputs:
        ----------
            gamma_air: float. piece-wise heat capacity ratio

        """

        # Check temperature range:
        if temp <= 600:
            return self.gamma_air_perfect
        elif 600 < temp <= 900:
            return self.gamma_air_comp
        elif 900 < temp <= 1400:
            return self.gamma_air_exp
        elif 1400 < temp <= 2000:
            return self.gamma_air_comb
        else:
            print('Another heat capacity ratio function should be chosen for temperatures above 2000K.')
            return self.gamma_air_comb
