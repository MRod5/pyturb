"""
IsentropicFlow class:
---------------------

Isentropic flow relations given an Ideal Gas.


MRodriguez 2020

"""

from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
import numpy as np


class IsentropicFlow(object):
    """
    IsentropicFlow:
    ---------------

    Isentropic flow relations for an Ideal Gas, both Perfect (constant cp, cv and
    gamma) or Semiperfect (cp, cv, gamma as a function of temperature).

    IsentropicFlow can inherit from PerfectIdealGas or SemiperfectIdealGas.

    Content:
    --------
        + sound_speed: Local speed of sound, assuming adiabatic flow
        + mach_number: Mach number for a local static temperature and velocity
        + stagnation_static_relation: Relation between stagnation temperature and
            static temperature as a function of Mach: Tt/T = 1+(gamma-1)/2*M**2
        + stag_temp_from_mach: Stagnation temperature given the stagnation to
            static relation and the static temperature: Tt = T * (Tt/T)
        + stag_temp_from_vel: Stagnation temperature given velocity of the flow
            Tt = T + v**2/cp/2
        + stag_pressure_from_mach: Stagnation pressure given the stagnation to
            static relation and the static pressure: pt = p * (Tt/T)**(gamma/(gamma-1))
        + impact_pressure_from_mach: Difference between the stagnation pressure and
            the static pressure.
        + dynamic_pressure: Dynamic pressure from static density and flow velocity.
            If M<0.3 incompressible flow may be consideed and pt = p + qi
        + dynamic_pressure_from_statpress: Dynamic pressure from static pressue and
            Mach number.
        + vel_from_stag_temp: Velocity of the flow given the stagnation and static
            temperatures
        + vel_from_mach: Veocity of the flow given the Mach number and the static temperature
        + stag_density_from_mach: Stagnation density given the stagnation to
            static relation and the static density: rhot = rho * (Tt/T)**(1/(gamma-1))


    """

    def __init__(self, fluid):
        if not(isinstance(fluid, PerfectIdealGas) or isinstance(fluid, SemiperfectIdealGas)):
            # Check the flow is a Perfect or a Semiperfect gas fom pyturb
            raise TypeError("Object must be PerfectIdealGas, SemiperfectIdealGas. Instead received {}".format(fluid))
        
        for attr in ['cp', 'gamma', 'Rg']:
            # Explicitly check for needed attributes
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
            # Avoid negative temperatures
            # TODO: Maybe allow unit change to international system
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
                # If Semiperfect gas a temperature must be provided
                raise ValueError("A valid static temperature must be provided")
            else:
                # Gamma for a Perfect gas
                gamma = self.fluid.gamma()
        else:
            # Gamma for semiperfect gas
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
            # Avoid negative temperatures
            # TODO: Maybe allow unit change to international system
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
            # Avoid negative temperatures
            # TODO: Maybe allow unit change to international system
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
                # If Semiperfect gas a temperature must be provided
                raise ValueError("A valid static temperature must be provided")
            else:
                # Gamma for a Perfect gas
                gamma = self.fluid.gamma()
                T = None
        else:
            # Gamma for a Semiperfect gas
            gamma = self.fluid.gamma(static_temperature)

            if not 0<static_temperature<6000:
                # Avoid negative temperatures
                # TODO: Maybe allow unit change to international system
                raise ValueError("Static temperature {}K out of bounds [0-6000]K".format(static_temperature))
            else:
                T = static_temperature

        if not 0<static_pressure:
            # Avoid negative pressure
            # TODO: Maybe allow unit change to international system
            raise ValueError("Static pressure {}Pa cannot be negative".format(static_pressure))
        else:
            p = static_pressure

        Tt_T = self.stagnation_static_rel(mach, T)
        pt = p * (Tt_T)**(gamma/(gamma-1))

        return pt


    def impact_pressure_from_mach(self, mach, static_pressure, static_temperature=None):
        """
        Impact pressure given the Mach number and the local static pressure.

        Inputs:
        -------
            mach: float. Mach number of the fluid [non-dimensional]
            static_pressure: float. Local static pressure [Pa]
            static_temperature: float. Static temperature of the gas [K]

        Outputs:
        --------
            qi: float. Impact pressure [Pa]

        """

        if not 0<static_pressure:
            # Avoid negative pressure
            raise ValueError("Static pressure {}Pa cannot be negative".format(static_pressure))
        else:
            p = static_pressure

        pt = self.stag_pressure_from_mach(mach, p, static_temperature)        
        qi = pt - p
        
        return qi


    def dynamic_pressure(self, vel, static_density):
        """
        Calculates the dynamic pressure of the flow, provided its static density
        and the velocity of the flow.
        
        Note that if the flow can be considered incompressible (isochoric), thus
        the Mach number of the flow is less than 0.3, then the stagnation pressure
        can be obtained as the static pressure plus the dynamic pressure:
            if M<0.3, then t = p + qi
        
        Inputs:
        -------
            vel: float. Flow velocity [m/s]
            static_density: float. Local static density [kg/m**3]

        Outputs:
        --------
            qi: float. Dynamic pressure of the flow [Pa]

        """

        if not 0<static_density:
            # Avoid negative density
            raise ValueError("Static density {}kg/m**3 cannot be negative".format(static_density))
        else:
            rho = static_density

        qi = 0.5*rho*vel**2

        return qi
        

    def dynamic_pressure_from_statpress(self, mach, static_pressure, static_temperature=None):
        """
        Calculates the dynamic pressure of the flow, provided its static pressure
        and the Mach number.
        
        Note that if the flow can be considered incompressible (isochoric), thus
        the Mach number of the flow is less than 0.3, then the stagnation pressure
        can be obtained as the static pressure plus the dynamic pressure:
            if M<0.3, then t = p + qi
        
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
            qi: float. Dynamic pressure of the flow [Pa]

        """

        if static_temperature is None:
            if not isinstance(self.fluid, PerfectIdealGas):
                # Gamma for a Perfect gas
                raise ValueError("A valid static temperature must be provided")
            else:
                # If Perfect gas a temperature must be provided
                gamma = self.fluid.gamma()
        else:
            # If Semiperfect gas a temperature must be provided
            gamma = self.fluid.gamma(static_temperature)

        if not 0<static_pressure:
            # Avoid negative pressure
            raise ValueError("Static pressure {}Pa cannot be negative".format(static_pressure))
        else:
            p = static_pressure

        qi = gamma/2 * p * mach**2

        return qi


    def vel_from_stag_temp(self, stagnation_temperature, static_temperature):
        """
        Flow velocity from the kinetic energy of the fluid, as the difference
        between the stagnation and static temperatures.

        Inputs:
        -------
            stagnation_temperature: float. Stagnation temperature of the gas [K]
            static_temperature: float. Static temperature of the gas [K]

        Outputs:
        --------
            vel: float. Flow velocity [m/s]

        """

        if not (0<static_temperature<6000 and 0<stagnation_temperature<6000):
            # Avoid negative temperatures
            raise ValueError("A valid temperatures must be provided.")
        elif not stagnation_temperature>static_temperature:
            # Avoid negative sqrt
            raise ValueError("Stagnation temperature must be greater than static temperature.")
        else:
            T = static_temperature
            Tt = stagnation_temperature
            

        difcpT = (self.fluid.cp(Tt)*Tt - self.fluid.cp(T)*T)
        vel = np.sqrt(2*difcpT)
        
        return vel


    def vel_from_mach(self, mach, static_temperature):
        """
        Flow velocity from the Mach number of the fluid.

        Inputs:
        -------
            mach: float. Mach number [dimensionless]
            static_temperature: float. Static temperature of the gas [K]


        Outputs:
        --------
            vel: float. Flow velocity [m/s]

        """

        if not 0<static_temperature<6000:
            # Avoid negative temperature
            raise ValueError("Static temperature {}K out of bounds [0-6000]K".format(static_temperature))
        else:
            T = static_temperature

        a = self.sound_speed(T)
        vel = a*mach

        return vel


    def stag_density_from_mach(self, mach, static_density, static_temperature=None):
        """
        Calculates the stagnation density given the Mach number and the local 
        static density.

        In case a semiperfect gas is used (gamma is a function of the temperature)
        a temperature must be provided to gather the heat capacity ratio.

        Inputs:
        -------
            mach: float. Mach number of the fluid [non-dimensional]
            static_density: float. Local static density [kg/m**3]
            static_temperature: float. Static temperature of the gas [K]
                If no static temperature is provided T_ref=298.15K is 
                taken if fluid is Perfect gas, otherwise a valid static 
                temperature must be provided.

        Outputs:
        --------
            rhot: float. stagnation density [Pa]

        """

        if static_temperature is None:
            if not isinstance(self.fluid, PerfectIdealGas):
                raise ValueError("A valid static temperature must be provided")
            else:
                # If Perfect gas a temperature must be provided
                gamma = self.fluid.gamma()
                T = None
        else:
            # If Semiperfect gas a temperature must be provided
            gamma = self.fluid.gamma(static_temperature)

            if not 0<static_temperature<6000:
                # Avoid negative temperatures
                raise ValueError("Static temperature {}K out of bounds [0-6000]K".format(static_temperature))
            else:
                T = static_temperature

        if not 0<static_density:
            # Avoid negative density
            raise ValueError("Static density {}kg/m**3 cannot be negative".format(static_density))
        else:
            rho = static_density

        Tt_T = self.stagnation_static_rel(mach, T)
        rhot = rho * (Tt_T)**(1/(gamma-1))

        return rhot

    