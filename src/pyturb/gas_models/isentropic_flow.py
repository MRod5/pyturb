"""
IsentropicFlow class:
---------------------

Isentropic flow relations


MRodriguez 2020

"""

from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
import numpy as np


class IsentropicFlow(object):
    """
    Isentropic flow relations.

    """

    def __init__(self, fluid):
        if not(isinstance(fluid, PerfectIdealGas) or isinstance(fluid, SemiperfectIdealGas)):
            raise TypeError("Object must be PerfectIdealGas, SemiperfectIdealGas. Instead received {}".format(fluid))
        
        for attr in ['cp', 'gamma']:
            if not hasattr(fluid, attr):
                raise AttributeError("Attribute {} not found in fluid object".format(attr))

        self.fluid = fluid
        return


    def sound_speed(self, static_temperature):
        """
        Local speed of sound, assuming adiabatic flow.
         
        Inputs:
        -------
            + static_temperature: float. Static temperature of the gas [K]
            
        Outputs:
        --------
            + a: speed of sound at a given static temperature, gamma and cp [m/s]
        """

        if not 0<static_temperature<6000:
            raise ValueError("Static temperature {}K out of bounds [0-6000]K".format(static_temperature))
        else:
            T = static_temperature

        a = np.sqrt(self.fluid.gamma(T) * self.fluid.Rg * T)

        return a


    def mach_number(self, velocity, static_temperature):
        """        
        Calculates the Mach number for a local static temperature condition and a given fluid velocity
        
        Inputs:
        -------
            velocity: float. Fluid velocity [m/s]
            static_temperature: float. Static temperature of the fluid [K]
        
        Outputs:
        --------
            mach: float. Mach number of the fluid [dimensionless]
            
        """
        mach = velocity / self.sound_speed(static_temperature)
        return mach


    def stagnation_static_rel(self, mach, static_temperature=None):
        """
        Calculates the stagnation to static relation given the Mach number of an
        isentropic flow.

        Inputs:
        -------
            mach: float. Mach number of the fluid [dimensionless]
            static_temperature: float. Static temperature of the gas [K]
                If no static temperature is provided T_ref=298.15K is 
                taken if fluid is Perfect gas, otherwise a valid static 
                temperature must be provided.

        Outputs:
        --------
            Tt_T: float. stagnation to static ratio [dimensionless]

        """

        if static_temperature is None:
            if not isinstance(self.fluid, PerfectIdealGas):
                raise ValueError("A valid static temperature must be provided")
            else:
                gamma = self.fluid.gamma()
        else:
            gamma = self.fluid.gamma(static_temperature)
        
        if 0>mach:
            raise ValueError("Wrong Mach number: {}".format(mach))

        Tt_T = 1 + (gamma - 1)/2 * mach**2
        return Tt_T
    
    
    def stag_temp_from_mach(self, mach, static_temperature):
        """
        Calculates the stagnation temperature of an isentropic flow given the
        Mach number and the static temperature of the flow.
        
        Inputs:
        -------
            mach: float. Mach number of the fluid [dimensionless]
            static_temperature: float. Static temperature of the gas [K]
                If no static temperature is provided T_ref=298.15K is 
                taken if fluid is Perfect gas, otherwise a valid static 
                temperature must be provided.

        Outputs:
        --------
            Tt: float. stagnation temperature [K]

        """

        if not 0<static_temperature<6000:
            raise ValueError("Static temperature {}K out of bounds [0-6000]K".format(static_temperature))
        else:
            T = static_temperature
        
        if 0>mach:
            raise ValueError("Wrong Mach number: {}".format(mach))

        Tt_T = self.stagnation_static_rel(mach, T)

        Tt = Tt_T * T
        return Tt


    def stag_temp_from_vel(self, velocity, static_temperature):
        """
        Calculates the stagnation temperature with the local velocity of the 
        fluid and the static temperature.

        Inputs:
        -------
            velocity: float. Fluid speed [m/s]
            static_temperature: float. Static temperature of the gas [K]
                If no static temperature is provided T_ref=298.15K is 
                taken if fluid is Perfect gas, otherwise a valid static 
                temperature must be provided.

        Outputs:
        --------
            Tt: float. stagnation temperature [K]
        """

        if not 0<static_temperature<6000:
            raise ValueError("Static temperature {}K out of bounds [0-6000]K".format(static_temperature))
        else:
            T = static_temperature

        Tt = T + velocity**2/self.fluid.cp(T)/2

        return Tt

    
    def stag_pressure_from_mach(self, mach, static_pressure, static_temperature=None):
        """
        Calculates the stagnation pressure given the Mach number and the local 
        static pressure.

        In case a semiperfect gas is used (gamma is a function of the temperature)
        a temperature must be provided to gather the heat capacity ratio.

        Inputs:
        -------
            mach: float. Mach number of the fluid [non-dimensional]
            static_pressure: float. Local static pressure [Pa]
            static_temperature: float. Static temperature of the gas [K]
                If no static temperature is provided T_ref=298.15K is 
                taken if fluid is Perfect gas, otherwise a valid static 
                temperature must be provided.

        Outputs:
        --------
            pt: float. stagnation pressure [Pa]

        """

        if static_temperature is None:
            if not isinstance(self.fluid, PerfectIdealGas):
                raise ValueError("A valid static temperature must be provided")
            else:
                gamma = self.fluid.gamma()
                T = None
        else:
            gamma = self.fluid.gamma(static_temperature)

            if not 0<static_temperature<6000:
                raise ValueError("Static temperature {}K out of bounds [0-6000]K".format(static_temperature))
            else:
                T = static_temperature

        if not 0<static_pressure:
            raise ValueError("Static pressure {}Pa cannot be negative".format(static_pressure))
        else:
            p = static_pressure

        Tt_T = self.stagnation_static_rel(mach, T)
        pt = p * (Tt_T)**(gamma/(gamma-1))

        return pt