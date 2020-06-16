"""
Nozzle:
-------

Generic nozzle control volume. Extends from ControlVolume.

Implements 30 different thermodynamic properties and variables of the 
control volume.


MRodriguez 2020

"""

from pyturb.power_plant.control_volume import ControlVolume
from pyturb.gas_models.isentropic_flow import IsentropicFlow
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas


class Nozzle(ControlVolume):
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


    ## Initializes the nozzle with the exit properties of a CV:
    def initialize_from_cv(self):
        """
        Initialize nozzle from the exit properties of last CV.
        """

        return


    ## Collects basic inputs of a nozzle:
    def initialize_nozzle(self, mflowe, pet, Tet, adiab_efficiency=1, Ae=None):
        """
        Set basic inputs of a generic nozzle:
            + 
            + 
        """

        self._mflow_e = mflowe
        self._p_et = pet
        self._T_et = Tet
        self._adiab_efficiency = adiab_efficiency
        self._A_e = Ae

        self.solve_basic_properties()

#        self.solve_adapted_nozzle()
#        self.solve_generic_nozzle()
#        self.solve_critical_convergent_nozzle()

        return None


    def solve_basic_properties(self):
        """
        """
        
        self._mflow_s = self.mflow_e
        self._delta_massflow = 0
        self._T_st = self.T_et
        self._h_et = self.fluid.cp(self.T_et) * self.T_et
        self._h_st = self.fluid.cp(self.T_st) * self.T_st
        self._q_se = 0
        self._w_se = 0

        if self.A_e is None:
            self._T_e = np.nan
            self._p_e = np.nan
            self._rho_e = np.nan
            self._rho_et = np.nan
            self._vel_e = np.nan
            self._h_e = np.nan
            self._ekin_e = np.nan
            self._mach_e = np.nan
        else:
            gamma_to = self.fluid.gamma(self.T_et)
            var_aux = self.mflow_e/self.A_e*np.sqrt(self.fluid.Rg/gamma_to) * np.sqrt(self.T_et)/self.p_et
            var_aux = var_aux ** (-2*(gamma_to-1)/(gamma_to+1))

            Me = symbols('Me')
            ec1 = Eq(var_aux, (Me+(1+(gamma_to-1)/2*Me**2) ) )

            M_e = solveset(ec1, Me, domain=S.Reals)
            M_e_ = list(M_e)
            M_e_.sort(reverse=True)

            mach_e_value = M_e_[0] if M_e_[0]>0 else None

            if mach_e_value is None:
                raise ValueError("Mach number at entrance of the nozzle is not possitive: {0}".format(mach_e_value))
            else:
                self._mach_e = float(mach_e_value)
            
            self._T_e = self.isent_flow.stat_temp_from_mach(self.mach_e, self.T_et)
            self._p_e = self.isent_flow.stat_pressure_from_mach(self.mach_e, self.p_et, self.T_et)
            self._rho_e = self.p_e / self.T_e / self.fluid.Rg
            self._rho_et = self.isent_flow.stag_density_from_mach(self.mach_e, self.rho_e, self.T_e)
            self._vel_e = self.isent_flow.vel_from_mach(self.mach_e, self.T_e)
            self._h_e = self.fluid.cp(self.T_e) * self.T_e
            self._ekin_e = 0.5 * self.vel_e**2

        return


    def solve_adapted_nozzle(self, ps, As=None):
        """
        """
        return
    

    def solve_critical_convergent_nozzle(self):
        """
        """

        return


    def solve_critical_condi_nozzle(self):
        """
        """

        return


    def solve_generic_nozzle(self):

        """
        """

        return

