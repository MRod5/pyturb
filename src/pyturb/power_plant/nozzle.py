"""
Nozzle:
-------

Generic nozzle control volume. Extends from ControlVolume.

Implements 30 different thermodynamic properties and variables of the 
control volume, plus critical conditions at the nozzle (A_star, T_star, p_star).


MRodriguez 2020

"""

from pyturb.power_plant.control_volume import ControlVolume
from pyturb.power_plant.intake import Intake
from pyturb.gas_models.isentropic_flow import IsentropicFlow
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
from sympy import symbols, Eq, solveset, S
import numpy as np
import warnings

class Nozzle(ControlVolume):
    """
    Nozzle:
    -------

    Defines a generic nozzle stage control volume, between one entrance and one
    exit to the control volume. the nozzle may be solved considering:
        + Adapted discharge
        + Critical, convergent nozzle
        + Critical, con-di (Laval) nozzle
        + Generic nozzle

    Contents:
    ---------
        + initialize_nozzle: initializes the control volume.
        + solve_basic_properties: partially solves basic thermodynamic properties of the nozzle
        + solve_adapted_nozzle: Completes the solution considering adapted discharge of the nozzle
        +
        +

    A total of 34 variables of a control volume are defined:
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
    
    + 4 properties at the nozzle throat, considering choked nozzle:
        + A_star, p_star, T_star

    """
    def __init__(self, fluid, stage=None):
        """
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


    ## Critical conditions:
    @property
    def exit_regime(self):
        """
        Flow regime at the exit of the nozzle.
        """
        return self._exit_regime


    @property
    def T_s_star(self):
        """
        Static temperature for choked nozzle throat.
        """
        return self._T_s_star


    @property
    def p_s_star(self):
        """
        Static temperature for choked nozzle throat.
        """
        return self._p_s_star


    @property
    def A_star(self):
        """
        Static temperature for choked nozzle throat.
        """
        return self._A_star


    ## Initializes the nozzle with the exit properties of a CV:
    def initialize_from_cv(self):
        """
        Initialize nozzle from the exit properties of last CV.
        """

        return


    ## Collects basic inputs of a nozzle:
    def initialize_nozzle(self, mflow_e, pet, Tet, Ae=None):
        """
        Set basic inputs of a generic nozzle:
            + mflow_e: float. Mass flow at the entrance [kg/s]
            + pet: float. stagnation pressure at the entrance [Pa]
            + Tet: float. stagnation tmeperature at the entrance [K]
            + Ae: float. Area at the entrance. [m**2]. If no area is provided
              related properties are set to np.nan in solve_basic_properties().
        """

        self._mflow_e = mflow_e
        self._p_et = pet
        self._T_et = Tet
        self._A_e = Ae

        self.solve_basic_properties()

        return None


    def solve_basic_properties(self):
        """
        Solves basic thermodyamic properties and CV variables at the entrance of
        the nozzle.
        """
        
        # Mass flow
        self._mflow_s = self.mflow_e # By definition
        self._delta_massflow = 0     # By definition
        
        # Enthalpies:
        self._T_st = self.T_et  # All nozzles isoenthalpic
        self._h_et = self.fluid.cp(self.T_et) * self.T_et
        self._h_st = self.fluid.cp(self.T_st) * self.T_st

        # By definition, no work/heat is done
        self._q_se = 0
        self._w_se = 0

        if self.A_e is None:
            # If the area at the entrance is not provided all area-related properties are set to nan
            self._T_e = np.nan
            self._p_e = np.nan
            self._rho_e = np.nan
            self._rho_et = np.nan
            self._vel_e = np.nan
            self._h_e = np.nan
            self._ekin_e = np.nan
            self._mach_e = np.nan

        else:
            # If the area is provided, the mach number is calculated
            gamma_to = self.fluid.gamma(self.T_et)

            # Solve mach at the entrance with variable size iterator:
            # Subsonic solution:
            var_aux = self.mflow_e/self.A_e*np.sqrt(self.fluid.Rg/gamma_to) * np.sqrt(self.T_et)/self.p_et
            expon = -(gamma_to+1)/(gamma_to-1)/2

            mach_func = lambda M: M*(1+(gamma_to-1)/2*M**2)**(expon) -var_aux
            mach_solution, _, _ = num_iters.variable_step_roots(x0=0, func=mach_func, dxmax=.2, verbosity=True)

            self._mach_e = mach_solution
            
            # Static properties at the entrance:
            self._T_e = self.isent_flow.stat_temp_from_mach(self.mach_e, self.T_et)
            self._p_e = self.isent_flow.stat_pressure_from_mach(self.mach_e, self.p_et, self.T_et)
            self._rho_e = self.p_e / self.T_e / self.fluid.Rg
            self._rho_et = self.isent_flow.stag_density_from_mach(self.mach_e, self.rho_e, self.T_e)
            self._vel_e = self.isent_flow.vel_from_mach(self.mach_e, self.T_e)
            self._h_e = self.fluid.cp(self.T_e) * self.T_e
            self._ekin_e = 0.5 * self.vel_e**2

        return


def solve_from_static_exit_pressure(self, ps, As=None, adiab_efficiency=1, nozzle_type='con-di'):
        """
        Solve nozzle with static discharge pressure known (exit section).
        The nozzle is assumed to be adiabatic.
            + If As is None, the discharge area is calculated with the provided
              adiabatic efficiency
            + If As is provided, the corresponding adiabatic efficiency is
              calculated. If the efficiency is impossible a warning is raised
        
        Inputs:
        -------
            ps: float. Static pressure at the exit section (stage #9). [Pa]
            As: float. Area of the exit section (stage #9) of a nozzle. If no area is 
                provided, it is calculated assuming the adiabatic efficiency set in
                initialize_nozzle. Otherwise the efficiency is recalculated. A warning
                is raised if the adiabatic efficiency is unfeasible.
            adiab_effiency: float. Adiabatic efficiency of the nozzle. If the exit area is
                not provided, isentropic nozzle is assumed. Otherwise it is calculated
            nozzle_type: string. May be 'con-di' for Laval/convergent-divergent nozzle or
                'convergent' for a convergent nozzle. By default 'con-di' is selected.

        """

        # Nozzle type
        if nozzle_type.lower()=='condi' or nozzle_type.lower()=='laval':
            nozzle_type = 'con-di'
        elif nozzle_type.lower()=='con' or nozzle_type.lower()=='conv':
            nozzle_type = 'convergent'
        elif not nozzle_type.lower() in ['con-di', 'convergent']:
            warnings.warn('Unknown nozzle type: {}. Nozzle will be set to con-di (default)'.format(nozzle_type))

        # Store static pressure
        self._p_s = ps

        gamma_to = self.fluid.gamma(self.T_et)

        if As is None:
            # Calculate discharge area assuming adapted nozzle with a given adiabatic efficiency
            self._adiab_efficiency = adiab_efficiency

            Ts_Tet =(1 + self.adiab_efficiency*((self.p_s/self.p_et)**((gamma_to-1)/gamma_to)- 1))
            self._T_s = self.T_et * Ts_Tet

            if self.T_s <= 2/(gamma_to + 1)*self.T_st:
                self._exit_regime = 'supersonic'
                self._T_s_star = self.isent_flow.stat_temp_from_mach(1, self.T_st)
                self._p_s_star = self.isent_flow.stat_pressure_from_mach(1, self.p_st, self.T_st)
                self._A_star = self.A_s*self.mach_s*((gamma_to+1)/2/(1+(gamma_to-1)/2*self.mach_s**2))**((gamma_to+1)/2/(gamma_to-1))

            else:
                self._exit_regime = 'subsonic'
                self._T_s_star = np.nan
                self._p_s_star = np.nan
                self._A_star = np.nan

            self._vel_s = self.isent_flow.vel_from_stag_temp(self.T_st, self.T_s)
            self._mach_s = self.isent_flow.mach_number(self.vel_s, self.T_s)

            self._p_st = self.isent_flow.stag_pressure_from_mach(self.mach_s, self.p_s, self.T_s)

            self._rho_s = self.p_s / self.fluid.Rg / self.T_s
            self._A_s = self.mflow_s / self.rho_s / self.vel_s

            self._rho_st = self.isent_flow.stag_density_from_mach(self.mach_s, self.rho_s, self.T_s)
            self._ekin_s = 0.5 * self.vel_s**2
            self._h_s = self.h_st - self.ekin_s

        else:
            # As is provided, check the adiabatic efficiency
            self._A_s = As

            # Solve static temperature from mass flow 2nd order equation:
            aux_var = 2*self.fluid.cp(self.T_st)*(self.p_s*self.A_s/self.fluid.Rg/self.mflow_s)**2
            T9 = (-aux_var + np.sqrt(aux_var**2 + 4*aux_var*self.T_st))/2
            self._T_s = T9

            if self.T_s <= 2/(gamma_to + 1)*self.T_st:
                self._exit_regime = 'supersonic'
                self._T_s_star = self.isent_flow.stat_temp_from_mach(1, self.T_st)
                self._p_s_star = self.isent_flow.stat_pressure_from_mach(1, self.p_st, self.T_st)
                self._A_star = self.A_s*self.mach_s*((gamma_to+1)/2/(1+(gamma_to-1)/2*self.mach_s**2))**((gamma_to+1)/2/(gamma_to-1))

            else:
                self._exit_regime = 'subsonic'
                self._T_s_star = np.nan
                self._p_s_star = np.nan
                self._A_star = np.nan

            self._vel_s = self.isent_flow.vel_from_stag_temp(self.T_st, self.T_s)
            self._mach_s = self.isent_flow.mach_number(self.vel_s, self.T_s)
            self._p_st = self.isent_flow.stag_pressure_from_mach(self.mach_s, self.p_s, self.T_s)
            self._rho_s = self.p_s / self.fluid.Rg / self.T_s

            self._rho_st = self.isent_flow.stag_density_from_mach(self.mach_s, self.rho_s, self.T_s)
            self._ekin_s = 0.5 * self.vel_s**2
            self._h_s = self.h_st - self.ekin_s

            # Recalculate adiabatic efficiency
            adiab_eff = (self.T_s/self.T_et-1)/((self.p_s/self.p_et)**((gamma_to-1)/gamma_to)-1)
            self._adiab_efficiency = adiab_eff
            if not 0<=adiab_eff<=1:
                warnings.warn('Unfeasible nozzle adiabatic efficiency (ad_eff={0}) for adapted nozzle (ps={1}) and fixed area (As={2})'.format(self.adiab_efficiency, self.p_s, self.A_s), UserWarning)
            

        if self.exit_regime=='supersonic':
            if nozzle_type=='convergent':
                # Recalculate solution assuming critical, convergent nozzles
                self.solve_critical_convergent_nozzle(self.ps, adiab_efficiency=self.adiab_efficiency)
            else:
                # Calculate nozzle regime. Solution is iterated:
                expon = (gamma_to+1)/(2*(gamma_to-1))
                area_rela = self.A_star/self.A_s

                # Mach function
                mach_func = lambda M: M*((gamma_to + 1)/2/( 1+(gamma_to-1)/2*M**2 ))**(expon) - area_rela

                # Subsonic solution (pressure at wich choking occurs)
                subsonic_mach, _, _ = num_iters.variable_step_roots(x0=0.5, func=mach_func, dxmax=.2, verbosity=True)

                # Supersonic solution (pressure at wich supersonic nozzle is adapted).
                supersonic_mach, _, _ = num_iters.variable_step_roots(x0=1.5, func=mach_func, dxmax=.2, verbosity=True)

                self._pchoke = self.isent_flow.stat_pressure_from_mach(subsonic_mach, self.p_st, self.T_st)
                self._padapt = self.isent_flow.stat_pressure_from_mach(supersonic_mach, self.p_st, self.T_st)

        return
    

    def solve_critical_convergent_nozzle(self, ps=None, As=None, adiab_efficiency=1):
        """
        Critical nozzle (choked) with convergent geometry (thus critical area is
        the discharge area A_s). Solves the nozzle assuming critical conditions and:
            + ps: If ps is provided, the discharge area and adiabatic efficiency are
                calculated
            + As: if As is provided, the corresponding adiabatic efficiency and static
              pressure are calculated. If the efficiency is impossible a warning is raised
            + adiab_efficiency: If no ps nor As are provided, the adiabatic efficiency is
                used to calculate the nozzle. By default the nozzle is assumed isentropic.

        Inputs:
        -------
            ps: float. Static discharge pressure [Pa]
            As: float. Discharge, critical area [m**2]
            adiab_efficiency: float. Adiabatic efficiency of the nozzle. By default is 1 (isentropic)
        
        """

        gamma_to = self.fluid.gamma(self.T_et)

        # Exit mach number and static temperature:
        self._mach_s = 1
        self._exit_regime = 'supersonic'
        self._T_s = self.isent_flow.stat_temp_from_mach(1, self.T_st)
        self._vel_s = self.isent_flow.sound_speed(self.T_s)

        if (As is None) and not (ps is None):
            # Store static pressure
            self._p_s = ps

            # Calculate the critical discharge area given the adiabatic efficiency
            self._p_st = self.isent_flow.stag_pressure_from_mach(self.mach_s, self.p_s, self.T_s)
            self._rho_s = self.p_s / self.fluid.Rg / self.T_s
            self._A_s = self.mflow_s / self.rho_s / self.vel_s

            self._rho_st = self.isent_flow.stag_density_from_mach(self.mach_s, self.rho_s, self.T_s)
            self._ekin_s = 0.5 * self.vel_s**2
            self._h_s = self.h_st - self.ekin_s

            # Recalculate adiabatic efficiency
            adiab_eff = (self.T_s/self.T_et-1)/((self.p_s/self.p_et)**((gamma_to-1)/gamma_to)-1)
            self._adiab_efficiency = adiab_eff
            if not 0<=adiab_eff<=1:
                warnings.warn('Unfeasible nozzle adiabatic efficiency (ad_eff={0}) for adapted nozzle (ps={1}) and fixed area (As={2})'.format(self.adiab_efficiency, self.p_s, self.T_s), UserWarning)
            
        elif not (As is None) and (ps is None):
           # As is provided
            self._A_s = As

            self._rho_s = self.mflow_s / self.A_s / self.vel_s
            self._rho_st = self.isent_flow.stag_density_from_mach(self.mach_s, self.rho_s, self.T_s)
            self._p_s = self.rho_s * self.fluid.Rg * self.T_s
            self._p_st = self.isent_flow.stag_pressure_from_mach(self.mach_s, self.p_s, self.T_s)

            self._rho_st = self.isent_flow.stag_density_from_mach(self.mach_s, self.rho_s, self.T_s)
            self._ekin_s = 0.5 * self.vel_s**2
            self._h_s = self.h_st - self.ekin_s

            # Recalculate adiabatic efficiency
            adiab_eff = (self.T_s/self.T_et-1)/((self.p_s/self.p_et)**((gamma_to-1)/gamma_to)-1)
            self._adiab_efficiency = adiab_eff
            if not 0<=adiab_eff<=1:
                warnings.warn('Unfeasible nozzle adiabatic efficiency (ad_eff={0}) for adapted nozzle (ps={1}) and fixed area (As={2})'.format(self.adiab_efficiency, self.p_s, self.T_s), UserWarning)

        elif (As is None) and (ps is None):
            # Calculate ps form adiabatic efficiency
            self._adiab_efficiency = adiab_efficiency
            
            # Static pressure
            self._p_s = self.p_et * (1 + (self.T_s/self.T_et -1)/self.adiab_efficiency)**(gamma_to/(gamma_to-1))
            self._rho_s = self.p_s / self.fluid.Rg / self.T_s
            self._A_s = self.mflow_s / self.rho_s / self.vel_s
            self._rho_st = self.isent_flow.stag_density_from_mach(self.mach_s, self.rho_s, self.T_s)            
            self._p_st = self.isent_flow.stag_pressure_from_mach(self.mach_s, self.p_s, self.T_s)

            self._ekin_s = 0.5 * self.vel_s**2
            self._h_s = self.h_st - self.ekin_s

        else:
            warnings.warn('If the nozzle is critical, the area and the static pressure at the discharge section cannot be set at the same time. The static pressure at the discharge section will be dismissed and recalculated.')
           # As is provided
            self._A_s = As

            self._rho_s = self.mflow_s / self.A_s / self.vel_s
            self._rho_st = self.isent_flow.stag_density_from_mach(self.mach_s, self.rho_s, self.T_s)
            self._p_s = self.rho_s * self.fluid.Rg * self.T_s
            self._p_st = self.isent_flow.stag_pressure_from_mach(self.mach_s, self.p_s, self.T_s)

            self._rho_st = self.isent_flow.stag_density_from_mach(self.mach_s, self.rho_s, self.T_s)
            self._ekin_s = 0.5 * self.vel_s**2
            self._h_s = self.h_st - self.ekin_s

            # Recalculate adiabatic efficiency
            adiab_eff = (self.T_s/self.T_et-1)/((self.p_s/self.p_et)**((gamma_to-1)/gamma_to)-1)
            self._adiab_efficiency = adiab_eff
            if not 0<=adiab_eff<=1:
                warnings.warn('Unfeasible nozzle adiabatic efficiency (ad_eff={0}) for adapted nozzle (ps={1}) and fixed area (As={2})'.format(self.adiab_efficiency, self.p_s, self.T_s), UserWarning)

        # Critical conditions:
        self._p_s_star = self.p_s
        self._T_s_star = self.T_s
        self._A_star = self.A_s

        if self.A_s>=self.A_e:
            warnings.warn('Exit area must be smaller than entry area in a convergent, critical nozzle. Results may be unfeasible')

        return


    def solve_generic_nozzle(self, ps=None, As=None, nozzle_type='con-di', adiab_efficiency=1):
        """
        Generic nozzle solver. Solves the nozzle assuming the discharge static
        pressure is not the ambient pressure (static)

        """

        # Nozzle type
        if nozzle_type.lower()=='condi' or nozzle_type.lower()=='laval':
            nozzle_type = 'con-di'
        elif nozzle_type.lower()=='con':
            nozzle_type = 'convergent'
        elif not nozzle_type.lower() in ['con-di', 'convergent']:
            warnings.warn('Unknown nozzle type: {}. Nozzle will be set to con-di (default)'.format(nozzle_type))
        

        gamma_to = self.fluid.gamma(self.T_et)

        if As is None:
            # Calculate discharge area assuming known adiabatic efficiency and static pressure
            # Store static pressure and efficiency:
            self._p_s = ps
            self._adiab_efficiency = adiab_efficiency

            Ts_Tet =(1 + self.adiab_efficiency*((self.p_s/self.p_et)**((gamma_to-1)/gamma_to)- 1))
            self._T_s = self.T_et * Ts_Tet

            if self.T_s <= 2/(gamma_to + 1)*self.T_st:
                self._exit_regime = 'supersonic'
            else:
                self._exit_regime = 'subsonic'

            self._vel_s = self.isent_flow.vel_from_stag_temp(self.T_st, self.T_s)
            self._mach_s = self.isent_flow.mach_number(self.vel_s, self.T_s)

            self._p_st = self.isent_flow.stag_pressure_from_mach(self.mach_s, self.p_s, self.T_s)

            self._rho_s = self.p_s / self.fluid.Rg / self.T_s
            self._A_s = self.mflow_s / self.rho_s / self.vel_s

            self._T_s_star = self.isent_flow.stat_temp_from_mach(1, self.T_st)
            self._p_s_star = self.isent_flow.stat_pressure_from_mach(1, self.p_st, self.T_st)
            self._A_star = self.A_s*self.mach_s*((gamma_to+1)/2/(1+(gamma_to-1)/2*self.mach_s**2))**((gamma_to+1)/2/(gamma_to-1))

            self._rho_st = self.isent_flow.stag_density_from_mach(self.mach_s, self.rho_s, self.T_s)
            self._ekin_s = 0.5 * self.vel_s**2
            self._h_s = self.h_st - self.ekin_s

        else:
            # XXX OJO CAMBIAR, AQUI SE CONOCE EL AREA Y LA EFICIENCIA PORQUE CONSTRUYES LA TOBERA, DESPEJAS PS
            # As is provided, check the adiabatic efficiency
            self._A_s = As

            # Solve static temperature from mass flow 2nd order equation:
            aux_var = 2*self.fluid.cp(self.T_st)*(self.p_s*self.A_s/self.fluid.Rg/self.mflow_s)**2
            T9 = (-aux_var + np.sqrt(aux_var**2 + 4*aux_var*self.T_st))/2
            self._T_s = T9

            if self.T_s <= 2/(gamma_to + 1)*self.T_st:
                self._exit_regime = 'supersonic'
            else:
                self._exit_regime = 'subsonic'

            self._vel_s = self.isent_flow.vel_from_stag_temp(self.T_st, self.T_s)
            self._mach_s = self.isent_flow.mach_number(self.vel_s, self.T_s)
            self._p_st = self.isent_flow.stag_pressure_from_mach(self.mach_s, self.p_s, self.T_s)
            self._rho_s = self.p_s / self.fluid.Rg / self.T_s

            self._T_s_star = self.isent_flow.stat_temp_from_mach(1, self.T_st)
            self._p_s_star = self.isent_flow.stat_pressure_from_mach(1, self.p_st, self.T_st)
            self._A_star = self.A_s*self.mach_s*((gamma_to+1)/2/(1+(gamma_to-1)/2*self.mach_s**2))**((gamma_to+1)/2/(gamma_to-1))

            self._rho_st = self.isent_flow.stag_density_from_mach(self.mach_s, self.rho_s, self.T_s)
            self._ekin_s = 0.5 * self.vel_s**2
            self._h_s = self.h_st - self.ekin_s

            # Recalculate adiabatic efficiency
            adiab_eff = (self.T_s/self.T_et-1)/((self.p_s/self.p_et)**((gamma_to-1)/gamma_to)-1)
            self._adiab_efficiency = adiab_eff
            if not 0<=adiab_eff<=1:
                warnings.warn('Unfeasible nozzle adiabatic efficiency (ad_eff={0}) for adapted nozzle (ps={1}) and fixed area (As={2})'.format(self.adiab_efficiency, self.p_s, self.T_s), UserWarning)
            

        if self.exit_regime=='supersonic' and nozzle_type=='convergent':
            # If the nozzle is convergent and supersonic, maximum mach number is 1
            # From the section where the nozzle is choked to the exit, the nozzle acts as a diffuser
            vel_star = self.isent_flow.vel_from_mach(1, self.T_s_star)
            
            supersonic_convergent_nozzle = Intake(self.fluid)
            supersonic_convergent_nozzle.initialize_intake(self.p_s_star, self.T_s_star,
            vel_star, self.A_star, adiab_efficiency=self.adiab_efficiency, As=self.A_s)


        return

