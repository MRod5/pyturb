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

class Intake(ControlVolume):
    """
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


    ## Collects basic inputs of an intake
    def initialize_intake(self, pe, Te, ve, Ae, adiab_efficiency=1, As=None):
        """
        Set basic inputs of a generic intake (diffuser):
            + Static pressure, static temperature and flow velocity at the entrance
            + Adiabatic efficiency of the diffuser
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
        """

        self._rho_e = self.p_e / self.T_e / self.fluid.Rg
        self._mach_e = self.isent_flow.mach_number(self.vel_e, self.T_e)
        self._T_et = self.isent_flow.stag_temp_from_mach(self.mach_e, self.T_e)
        self._p_et = self.isent_flow.stag_pressure_from_mach(self.mach_e, self.p_e, self.T_e)
        self._rho_et = self.isent_flow.stag_density_from_mach(self.mach_e, self.rho_e, self.T_e)

        self._mflow_e = self.rho_e * self.A_e * self.vel_e
        self._mflow_s = self.mflow_e
        self._delta_massflow = self.mflow_s - self.mflow_e

        self._h_e = self.fluid.cp(self.T_e) * self.T_e
        self._ekin_e = 0.5 * self.vel_e ** 2
        self._h_et = self.h_e + self.ekin_e
        self._h_st = self.h_et

        self._T_st = self.T_et

