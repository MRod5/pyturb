#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
isentropic_gas:
---------------

This file contains isentropic flow relations (stagnation temperatures and pressures), for a given model of air 
(thermally perfect of thermally real, see air_model.py).
    
Content:
--------
    + IsentropicGas: isentropic relations
              
Part of pyturb. Tests pending.

April 2017

MRodriguez
"""

import numpy as np
from air_model import RealGas, IdealGas


class IsentropicGas:
    """
    Isentropic flow relations for a given model of air (thermally perfect of thermally real).
    
    Content:
    --------
    + sound_speed: local sound speed
    + mach_number: Mach number of the fluid
    + stagnation_static_temperature_relation: relation between the stagnation and static temperature in terms of the 
                                              Mach number
    + stagnation_temperature_from_vel: stagnation temperature obtained with static temperature and fluid velocity
    + stagnation_temperature_from_mach: stagnation temperature obtained with static temperature and Mach number
    + stagnation_static_pressure_relation: relation between the stagnation and static pressure in terms of the
                                           Mach number
    + stagnation_pressure_from_mach: stagnation pressure obtained with static pressure and Mach number
    + velocity_from_stagnation_temperature: flow velocity from stagnation and static pressure
    + mach_from_stagnation_temperature: Mach number of the gas from stagnation and static temperature
    + mach_from_stagnation_pressure: Mach number of the gas from stagnation and static pressure
    + stagnation_enthalpy: stagnation specific enthalpy from static enthalpy en kinetic specific energy
    + kinetic_energy_from_enthalpy: specific kinetic energy of the fluid from the stagnation and static enthalpy
    + velocity_from_enthalpy: flow velocity from stagnation and static specific enthalpy
    
        
    """

    def __init__(self, selected_cp_air_model='ideal', selected_gamma_air_model='ideal'):
        """
        Initializes IsentropicGas. Switches between RealGas or IdealGas from air_model.py depending on the inputs:
        + Inputs:
        ---------
            selected_cp_air_model: 'ideal' or any real model defined at air_model.RealGas
            selected_gamma_air_model: 'ideal' or any real model defined at air_model.RealGas
         
        """
        if type(selected_cp_air_model) == str:
            if selected_cp_air_model.lower() == 'ideal':
                self.selected_cp_air_model = 'ideal'
                if selected_gamma_air_model != 'ideal':
                    print('Gamma_air model switched to ideal')
                self.selected_cp_air_model = 'ideal'

                self.air_data = IdealGas()

            else:
                am = RealGas()
                available_cp_functions = am.cp_air_functions.keys()
                if selected_cp_air_model.lower() in available_cp_functions:
                    self.selected_cp_air_model = selected_cp_air_model.lower()

                    available_gamma_functions = am.gamma_air_functions.keys()
                    if selected_gamma_air_model.lower() in available_gamma_functions:
                        self.selected_gamma_air_model = selected_gamma_air_model.lower()
                        self.air_data = RealGas(cp_option=self.selected_cp_air_model,
                                                gamma_option=selected_gamma_air_model)
                    else:
                        print('Gamma_air model switched to default')
                        self.air_data = RealGas(cp_option=self.selected_cp_air_model)
                else:
                    print('Unknown selected cp real gas model. Switching to Ideal Gas')
                    self.air_data = IdealGas()
        else:
            print('Unknown selected cp model. Switching to ideal gas...')

            self.selected_cp_air_model = 'ideal'
            self.selected_gamma_air_model = 'ideal'
            self.air_data = IdealGas()

        # Specific gas constant:
        self.R_gas = self.air_data.R_air

    def _gamma_gas(self, static_temperature):
        """
        Gets the specific heat ratio depending on the selected air model:
        
        + Inputs:
        ---------
            static_temperature: float. Static temperature of the gas [K]
        
        + Outputs:
        ----------
            gamma_gas: float. Specific heat ratio [non-dimensional]
        """

        return self.air_data.gamma_air(static_temperature)

    def _cp_gas(self, static_temperature):
        """
        Gets the specific heat at constant pressure of the gas depending on the selected air model.

        + Inputs:
        ---------
            static_temperature: float. Static temperature of the gas [K]

        + Outputs:
        ----------
            cp_gas: float. Specific heat at constant pressure [J/kg/K]
        """

        return self.air_data.cp_air(static_temperature)

    def sound_speed(self, static_temperature):
        """
        Calculates local speed of sound.
         
        Inputs:
        -------
            + static_temperature: float. Static temperature of the gas [K]
            
        Outputs:
        --------
            + a: speed of sound at a given static temperature, gamma and cp [m/s]
         
        """
        return np.sqrt(self._gamma_gas(static_temperature) * self.R_gas * static_temperature)

    def mach_number(self, velocity, static_temperature):
        """        
        Calculates the Mach number for a local static temperature condition and a given fluid velocity
        
        Inputs:
        -------
            velocity: float. Fluid speed [m/s]
            static_temperature: float. Static temperature of the gas [K]
        
        Outputs:
        --------
            mach: float. Mach number of the fluid [non-dimensional]
            
        """
        mach = velocity / self.sound_speed(static_temperature)
        return mach

    def stagnation_static_temperature_relation(self, mach_number, static_temperature):
        """        
        Calculates the stagnation temperature to static temperature ratio by means of the mach number

        Inputs:
        -------
            mach: float. Mach number of the fluid [non-dimensional]
            static_temperature: float. Static temperature of the gas [K]

        Outputs:
        --------
            stag_stat_temp: float. stagnation temperature to static temperature ratio [non-dimensional]

        """

        return 1 + (self._gamma_gas(static_temperature) - 1) / 2 * mach_number ** 2

    def stagnation_temperature_from_vel(self, velocity, static_temperature):
        """
        Calculates the stagnation temperature with the local velocity of the fluid and the static temperature

        Inputs:
        -------
            velocity: float. Fluid speed [m/s]
            static_temperature: float. Static temperature of the gas [K]

        Outputs:
        --------
            temp_t: float. stagnation temperature [K]
        """

        return static_temperature + 0.5*velocity**2/self._cp_gas(static_temperature)

    def stagnation_temperature_from_mach(self, mach, static_temperature):
        """
        Calculates the stagnation temperature as a function of the local Mach number and the static temperature

        Inputs:
        -------
            mach: float. Mach number of the fluid [non-dimensional]
            static_temperature: float. Static temperature of the gas [K]

        Outputs:
        --------
            temp_t: float. stagnation temperature [K]
        """

        return static_temperature*self.stagnation_static_temperature_relation(mach, static_temperature)

    def stagnation_static_pressure_relation(self, mach, static_temperature):
        """
        Calculates the stagnation to static pressure relation with the local Mach number

        Inputs:
        -------
            mach: float. Mach number of the fluid [non-dimensional]
            static_pressure: float. Static pressure of the gas [Pa]

        Outputs:
        --------
            stag_stat_press_rel: float. stagnation to static pressure relation [non-dimensional]
        """

        return (self.stagnation_static_temperature_relation(mach, static_temperature))**(
            self._gamma_gas(static_temperature)/(self._gamma_gas(static_temperature) - 1))

    def stagnation_pressure_from_mach(self, mach, static_pressure, static_temperature):
        """
        Calculates the stagnation pressure as a function of the local Mach number and the static pressure.

        Inputs:
        -------
            mach: float. Mach number of the fluid [non-dimensional]
            static_pressure: float. Static pressure of the gas [Pa]

        Outputs:
        --------
            press_t: float. stagnation pressure [Pa]
        """

        return static_pressure * (self.stagnation_static_pressure_relation(mach, static_temperature))

    def velocity_from_stagnation_temperature(self, stagnation_temperature, static_temperature):
        """
        Calculates the velocity of the fluid provided its stagnation and static temperatures.
        
        + Inputs:
        ---------
            stagnation_temperature: float. Stagnation temperature of the gas [K]
            static_temperature: float. Static temperature of the gas [K]
        
        + Outputs:
        ----------
            velocity: float. Velocity of the fluid [m/s]
        """

        return np.sqrt(2*self._cp_gas(static_temperature)*(stagnation_temperature - static_temperature))

    def mach_from_stagnation_temperature(self, stagnation_temperature, static_temperature):
        """
        Calculates the Mach number provided the stagnation and static temperatures of the gas

        + Inputs:
        ---------
            stagnation_temperature: float. Stagnation temperature of the gas [K]
            static_temperature: float. Static temperature of the gas [K]

        + Outputs:
        ----------
            mach: float. Mach number of the fluid [non-dimensional]
        """

        mach = np.sqrt(2/(self._gamma_gas(static_temperature) - 1)*(stagnation_temperature / static_temperature - 1))
        return mach

    def mach_from_stagnation_pressure(self, stagnation_pressure, static_pressure, static_temperature):
        """
        Calculates the Mach number provided the stagnation and static pressures of the gas and the static temperature

        + Inputs:
        ---------
            stagnation_pressure: float. Stagnation pressure of the gas [Pa]
            static_pressure: float. Static temperature of the gas [Pa]
            static_temperature: float. Static temperature of the gas [K]

        + Outputs:
        ----------
            mach: float. Mach number of the fluid [non-dimensional]
        """
        gamma_g = self._gamma_gas(static_temperature)
        mach = np.sqrt(2 / (gamma_g - 1) * ((stagnation_pressure / static_pressure)**((gamma_g - 1)/gamma_g) - 1))
        return mach

    def stagnation_enthalpy(self, velocity, static_temperature):
        """
        Calculates the stagnation specific enthalpy from the static enthalpy and the specific kinetic energy

        Inputs:
        -------
            velocity: float. Fluid speed [m/s]
            static_temperature: float. Static temperature of the gas [K]

        Outputs:
        --------
            stagnation_enthalpy: float. stagnation specific enthalpy [J/kg]
        """

        return self._cp_gas(static_temperature)*static_temperature + 0.5*velocity**2

    @staticmethod
    def kinetic_energy_from_enthalpy(stagnation_temperature, static_temperature):
        """
        Calculates the velocity of the fluid with the stagnation and static enthalpy.

        Inputs:
        -------
            stagnation_temperature: float. Stagnation enthalpy of the gas [J/kg]
            static_temperature: float. Static enthalpy of the gas [J/kg]

        Outputs:
        --------
            kinetic_energy: float. kinetic specific energy [J/kg]
        """

        return stagnation_temperature - static_temperature

    def velocity_from_enthalpy(self, stagnation_temperature, static_temperature):
        """
        Calculates the velocity of the fluid with the stagnation and static enthalpy

        Inputs:
        -------
            stagnation_temperature: float. Stagnation enthalpy of the gas [J/kg]
            static_temperature: float. Static enthalpy of the gas [J/kg]

        Outputs:
        --------
            stagnation_enthalpy: float. stagnation specific enthalpy [J/kg]
        """

        return np.sqrt(2*self.kinetic_energy_from_enthalpy(stagnation_temperature, static_temperature))
