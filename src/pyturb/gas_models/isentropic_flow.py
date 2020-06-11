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