"""
Control Volume:
---------------

Class for a generic control volume, with static thermodynamic properties, stagnation
properties, Energy Conservation (First Law), adiabatic and polytropic efficiency (Second
Law) and Conservation of Mass (Continuity) for a generic control volume (CV).

A total of 30 variables of a control volume are defined:
    + 12 thermodynamic properties at the entrance of the CV:
    + 12 thermodynamic properties at the exit of the CV
    + 6 variables of the CV

MRodriguez 2020

"""

from abc import abstractmethod
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas



class ControlVolume(object):
    """
    Control Volume:
    ---------------

    A total of 30 variables of a control volume are defined:

    + 12 thermodynamic properties at the entrance of the CV:
        + Static properties:
            + p_e, T_e, rho_e, h_e (enthalpy)
        + Stagnation properties
            + p_et, T_et, rho_et, h_et (enthalpy)
        + Other:
            + vel_e (flow velocity), mach_e, ec_e (kin. energy), mflow_e (mass flow)

    + 12 thermodynamic properties at the exit of the CV
        + Static properties:
            + p_s, T_s, rho_s, h_s (enthalpy)
        + Stagnation properties
            + p_st, T_st, rho_st, h_st (enthalpy)
        + Other:
            + vel_s (flow velocity), mach_s, ec_s (kin. energy), mflow_s (mass flow)


    + 6 variables of the CV:
        + specific power:
            + q_se (heating/cooling power along the CV), w_se (mechanical power along the CV)
        + adiab_efficiency, A_e (entrance area), A_s (exit area), delta_massflow (mass flow variation along CV)

    A pyturb fluid class is needed to initialize a ControlVolume.



    """
    def __init__(self, fluid):
        if not(isinstance(fluid, PerfectIdealGas) or isinstance(fluid, SemiperfectIdealGas)):
            raise TypeError("PerfectIdealGas, SemiperfectIdealGas")
        self.fluid = fluid
        return

    # **************************
    # *** Entrance of the CV *** 
    # **************************
    # Static thermodynamic properties:
    @property
    @abstractmethod
    def p_e(self):
        """
        Static pressure at the entrance of the CV. [Pa]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def T_e(self):
        """
        Static temperature at the entrance of the CV. [K]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def rho_e(self):
        """
        Static density at the entrance of the CV. [kg/m**3]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def h_e(self):
        """
        Specific static enthalpy at the entrance of the CV. [W/kg/s]
        """
        raise NotImplementedError


    # Stagnation thermodynamic properties
    @property
    @abstractmethod
    def p_et(self):
        """
        Stagnation pressure at the entrance of the CV. [Pa]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def T_et(self):
        """
        Stagnation temperature at the entrance of the CV. [K]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def rho_et(self):
        """
        Stagnation density at the entrance of the CV. [kg/m**3]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def h_et(self):
        """
        Specific stagnation enthalpy at the entrance of the CV. [W/kg/s]
        """
        raise NotImplementedError


    # Velocity and kinetic energy:
    @property
    @abstractmethod
    def vel_e(self):
        """
        Flow velocity at the entrance of the CV. [m/s]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def mach_e(self):
        """
        Mach number at the entrance of the CV. [dimensionless]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def ec_e(self):
        """
        Specific kinetic energy at the entrance of the CV. [W/kg/s]
        """
        raise NotImplementedError
    
    
    # Mass flow and area:
    @property
    @abstractmethod
    def mflow_e(self):
        """
        Mass flow at the entrance of the CV. [kg/s]
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def A_e(self):
        """
        Area at the entrance of the CV. [m**2]
        """
        raise NotImplementedError


    # **********************
    # *** Exit of the CV *** 
    # **********************
        # Static thermodynamic properties:
    @property
    @abstractmethod
    def p_s(self):
        """
        Static pressure at the entrance of the CV. [Pa]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def T_s(self):
        """
        Static temperature at the entrance of the CV. [K]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def rho_s(self):
        """
        Static density at the entrance of the CV. [mg/m**3]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def h_s(self):
        """
        Specific static enthalpy at the entrance of the CV. [W/kg/s]
        """
        raise NotImplementedError


    # Stagnation thermodynamic properties
    @property
    @abstractmethod
    def p_st(self):
        """
        Stagnation pressure at the entrance of the CV. [Pa]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def T_st(self):
        """
        Stagnation temperature at the entrance of the CV. [K]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def rho_st(self):
        """
        Stagnation density at the entrance of the CV. [kg/m**3]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def h_st(self):
        """
        Specific stagnation enthalpy at the entrance of the CV. [W/kg/s]
        """
        raise NotImplementedError


    # Velocity and kinetic energy:
    @property
    @abstractmethod
    def vel_s(self):
        """
        Flow velocity at the entrance of the CV. [m/s]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def mach_s(self):
        """
        Mach number at the entrance of the CV. [dimensionless]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def ec_s(self):
        """
        Specific kinetic energy at the entrance of the CV. [W/kg/s]
        """
        raise NotImplementedError
    
    
    # Mass flow and area:
    @property
    @abstractmethod
    def mflow_s(self):
        """
        Mass flow at the entrance of the CV. [kg/s]
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def A_s(self):
        """
        Area at the entrance of the CV. [m**2]
        """
        raise NotImplementedError


    # ***************************
    # *** Transfer variables: *** 
    # ***************************
    @property
    @abstractmethod
    def w_se(self):
        """
        Speficic mechanical power along the CV. [W/mg/s]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def q_se(self):
        """
        Speficic heating power along the CV. [W/mg/s]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def delta_massflow(self):
        """
        Mass flow variation along the CV. [W/mg/s]
        """
        raise NotImplementedError


    @property
    @abstractmethod
    def adiab_efficiency(self):
        """
        Adiabatic efficiency along the CV. [dimensionless]
        """
        raise NotImplementedError