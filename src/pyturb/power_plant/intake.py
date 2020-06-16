"""
Intake:
-------

Generic intake (diffuser) control volume. Extends from ControlVolume.

Implements 30 different thermodynamic properties and variables of the 
control volume.


MRodriguez 2020

"""

from pyturb.power_plant.control_volume import ControlVolume
from pyturb.gas_models.isentropic_flow import IsentropicFlow
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
from sympy import symbols, Eq, solveset, S
import numpy as np


class Intake(ControlVolume):
    """
    Intake:
    -------

    Defines a generic intake stage control volume, between one entrance and one
    exit to the control volume.

    Contents:
    ---------
        + initialize_intake: initializes the control volume with the minimum number
            of thermodynamic properties and variables needed to solve the intake.
        + solve_intake: solves all thermodynamic properties and variables involved
            in the control volume.

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

    """
    def __init__(self, fluid, stage=None):
        """
        Initializes Intake and ControlVolume.
        """
        # Stage of the power plant:
        self.stage = stage
        
        # Initialize superclass methods
        super().__init__(fluid)

        # Isentropic flow:
        if not(isinstance(fluid, PerfectIdealGas) or isinstance(fluid, SemiperfectIdealGas)):
            # Check the flow is a Perfect or a Semiperfect gas fom pyturb
            raise TypeError("Object must be PerfectIdealGas, SemiperfectIdealGas. Instead received {}".format(fluid))
    
        self.isent_flow = IsentropicFlow(self.fluid)

            
        return None


    # **************************
    # *** Entrance of the CV *** 
    # **************************
    # Static thermodynamic properties:
    @property
    def p_e(self):
        """
        Static pressure at the entrance of the CV. [Pa]
        """
        return self._p_e
    

    @property
    def T_e(self):
        """
        Static temperature at the entrance of the CV. [K]
        """
        return self._T_e
    

    @property
    def rho_e(self):
        """
        Static density at the entrance of the CV. [kg/m**3]
        """
        return self._rho_e


    @property
    def h_e(self):
        """
        Specific static enthalpy at the entrance of the CV. [W/kg/s]
        """
        return self._h_e


    # Stagnation thermodynamic properties
    @property
    def p_et(self):
        """
        Stagnation pressure at the entrance of the CV. [Pa]
        """
        return self._p_et


    @property
    def T_et(self):
        """
        Stagnation temperature at the entrance of the CV. [K]
        """
        return self._T_et


    @property
    def rho_et(self):
        """
        Stagnation density at the entrance of the CV. [kg/m**3]
        """
        return self._rho_et


    @property
    def h_et(self):
        """
        Specific stagnation enthalpy at the entrance of the CV. [W/kg/s]
        """
        return self._h_et


    # Velocity and kinetic energy:
    @property
    def vel_e(self):
        """
        Flow velocity at the entrance of the CV. [m/s]
        """
        return self._vel_e


    @property
    def mach_e(self):
        """
        """
        return self._mach_e
    

    @property
    def ekin_e(self):
        """
        Specific kinetic energy at the entrance of the CV. [W/kg/s]
        """
        return self._ekin_e


    # Mass flow and area:
    @property
    def mflow_e(self):
        """
        Mass flow at the entrance of the CV. [kg/s]
        """
        return self._mflow_e


    @property
    def A_e(self):
        """
        Mass flow at the entrance of the CV. [kg/s]
        """
        return self._A_e


    # **********************
    # *** Exit of the CV *** 
    # **********************
    # Static thermodynamic properties:
    @property
    def p_s(self):
        """
        Static pressure at the entrance of the CV. [Pa]
        """
        return self._p_s


    @property
    def T_s(self):
        """
        Static temperature at the entrance of the CV. [K]
        """
        return self._T_s


    @property
    def rho_s(self):
        """
        Static density at the entrance of the CV. [mg/m**3]
        """
        return self._rho_s


    @property
    def h_s(self):
        """
        Specific static enthalpy at the entrance of the CV. [W/kg/s]
        """
        return self._h_s


    # Stagnation thermodynamic properties
    @property
    def p_st(self):
        """
        Stagnation pressure at the entrance of the CV. [Pa]
        """
        return self._p_st


    @property
    def T_st(self):
        """
        Stagnation temperature at the entrance of the CV. [K]
        """
        return self._T_st


    @property
    def rho_st(self):
        """
        Stagnation density at the entrance of the CV. [kg/m**3]
        """
        return self._rho_st


    @property
    def h_st(self):
        """
        Specific stagnation enthalpy at the entrance of the CV. [W/kg/s]
        """
        return self._h_st


    # Velocity and kinetic energy:
    @property
    def vel_s(self):
        """
        Flow velocity at the entrance of the CV. [m/s]
        """
        return self._vel_s


    @property
    def mach_s(self):
        """
        Mach number at the entrance of the CV. [dimensionless]
        """
        return self._mach_s


    @property
    def ekin_s(self):
        """
        Specific kinetic energy at the entrance of the CV. [W/kg/s]
        """
        return self._ekin_s
    
    
    # Mass flow and area:
    @property
    def mflow_s(self):
        """
        Mass flow at the entrance of the CV. [kg/s]
        """
        return self._mflow_s


    @property
    def A_s(self):
        """
        Area at the entrance of the CV. [m**2]
        """
        return self._A_s


    # ***************************
    # *** Transfer variables: *** 
    # ***************************
    @property
    def w_se(self):
        """
        Speficic mechanical power along the CV. [W/mg/s]
        """
        return self._w_se


    @property
    def q_se(self):
        """
        Speficic heating power along the CV. [W/mg/s]
        """
        return self._q_se


    @property
    def delta_massflow(self):
        """
        Mass flow variation along the CV. [W/mg/s]
        """
        return self._delta_massflow


    @property
    def adiab_efficiency(self):
        """
        Adiabatic efficiency along the CV. [dimensionless]
        """
        return self._adiab_efficiency

    
    ## Collects basic inputs of an intake
    def initialize_intake(self, pe, Te, ve, Ae, adiab_efficiency=1, As=None):
        """
        Sets basic inputs of a generic intake (diffuser) control volume:
            + Static pressure, static temperature and flow velocity at the entrance
            + Adiabatic efficiency of the diffuser (by default the intake is considered
              to be isentropic)
            + Exit area of the intake (if no exit area is provided, properties dependent
            on the area are set to np.nan).
        """

        self._A_e = Ae
        self._p_e = pe
        self._T_e = Te
        self._vel_e = ve
        self._adiab_efficiency = adiab_efficiency
        self._A_s = As

        self.solve_intake()

        return None


    ## Solves the intake at the entrance:
    def solve_intake(self):
        """
        Solves the thermodynamic properties and variables of the intake.

        Note that the followwing variables and properties are solved only if an exit area
        is provided in initialize_intake(), otherwise those variables are set to np.nan:
            vel_s, p_s, T_s, rho_s, ekin_s, h_s, mach_s
        
        If the exit area is provided, Sympy is used to solve the mach at the exit
        of the intake.
        """

        # Solve basic thermodyamic properties at the entrance:
        self._rho_e = self.p_e / self.T_e / self.fluid.Rg
        self._mach_e = self.isent_flow.mach_number(self.vel_e, self.T_e)
        self._T_et = self.isent_flow.stag_temp_from_mach(self.mach_e, self.T_e)
        self._p_et = self.isent_flow.stag_pressure_from_mach(self.mach_e, self.p_e, self.T_e)
        self._rho_et = self.isent_flow.stag_density_from_mach(self.mach_e, self.rho_e, self.T_e)

        # Mass flow along the control volume
        self._mflow_e = self.rho_e * self.A_e * self.vel_e
        # TODO: If the intake is used with a fan, bypass ratio may be obtained here
        self._mflow_s = self.mflow_e
        self._delta_massflow = self.mflow_s - self.mflow_e

        # Enthalpies:
        self._h_e = self.fluid.cp(self.T_e) * self.T_e
        self._ekin_e = 0.5 * self.vel_e ** 2
        self._h_et = self.h_e + self.ekin_e
        self._h_st = self.h_et

        # Thermodynamic properties at the exit fo the CV
        self._T_st = self.T_et

        gamma_d = self.fluid.gamma(self.T_et)
        self._p_st = self.p_e * ((self.T_st / self.T_e - 1)*self.adiab_efficiency + 1)**(gamma_d/(gamma_d-1))
        self._rho_st = self.rho_e * ((self.T_st / self.T_e - 1)*self.adiab_efficiency + 1)**(1/(gamma_d-1))

        # By definition, no work/heat is applied in a diffuser cv
        self._w_se = 0
        self._q_se = 0

        if self.A_s is None:
            # Is exit area is not provided:
            self._vel_s = np.nan
            self._p_s = np.nan
            self._T_s = np.nan
            self._rho_s = np.nan
            self._ekin_s = np.nan
            self._h_s = np.nan
            self._mach_s = np.nan
        else:
            # With the exit area, solve the mach number at the exit
            var_aux = self.mflow_s/self.A_s*np.sqrt(self.fluid.Rg/gamma_d) * np.sqrt(self.T_st)/self.p_st
            var_aux = var_aux ** (-2*(gamma_d-1)/(gamma_d+1))

            Ms = symbols('Ms')

            ec1 = Eq(var_aux, (Ms+(1+(gamma_d-1)/2*Ms**2) ) )

            # Dismiss negative solutions
            M_s = solveset(ec1, Ms, domain=S.Reals)
            M_s_ = list(M_s)
            M_s_.sort(reverse=True)
            mach_s_value = M_s_[0] if M_s_[0]>0 else None

            if mach_s_value is None:
                raise ValueError("Mach number at exit of the CV is not possitive: {0}".format(mach_s_value))
            else:
                self._mach_s = float(mach_s_value)
            
            # Calculate static properties at the exit of the CV
            self._T_s = self.isent_flow.stat_temp_from_mach(self.mach_s, self.T_st)
            self._p_s = self.isent_flow.stat_pressure_from_mach(self.mach_s, self.p_st, self.T_st)
            self._vel_s = self.isent_flow.vel_from_mach(self.mach_s, self.T_s)
            self._rho_s = self.p_s / self.T_s / self.fluid.Rg
            self._ekin_s = 0.5 * self.vel_s**2
            self._h_s = self.fluid.cp(self.T_s) * self.T_s

        return

