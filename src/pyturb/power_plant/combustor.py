"""
Combustor:
----------

Generic combustion chamber control volume. Extends from ControlVolume.

Implements 30 different thermodynamic properties and variables of the 
control volume.


MRodriguez 2020

"""

from pyturb.power_plant.control_volume import ControlVolume
from pyturb.gas_models.isentropic_flow import IsentropicFlow
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
import numpy as np
import warnings


class Combustor(ControlVolume):
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


    ## Combustion parameters
    @property
    def combustion_efficiency(self):
        """
        Heat of combustion absorption efficiency. [dimensionless]
        """
        return self._adiab_efficiency


    ## Initializes the combustion chamber with the exit properties of a CV:
    def initialize_from_cv(self, cv):
        """
        Initialize nozzle from the exit properties of last CV.
        """

        # Isentropic flow:
        if not(isinstance(cv, ControlVolume)):
            # The control volume from which the nozzle has to be initialized is not a ControlVolume object
            raise TypeError("Object must be of type ControlVolume, SemiperfectIdealGas. Instead received {}".format(cv))
    

        self._mflow_e = cv.mflow_s
        self._p_et = cv.p_st
        self._T_et = cv.T_st
        self._A_e = cv.A_s

        self.solve_basic_properties()

        return


    ## Collects basic inputs of a combustion chamber:
    def initialize_nozzle(self, mflow_e, pet, Tet, Ae=None, c=None, f=None, phi=None):
        """
        Set basic inputs of a generic nozzle:
            + mflow_e: float. Mass flow at the entrance [kg/s]
            + pet: float. stagnation pressure at the entrance [Pa]
            + Tet: float. stagnation temperature at the entrance [K]
            + Ae: float. Area at the entrance. [m**2]. If no area is provided
              related properties are set to np.nan in solve_basic_properties()
        """

        self._mflow_e = mflow_e
        self._p_et = pet
        self._T_et = Tet
        self._A_e = Ae
        self._c = c
        self._far = f
        self._equiv_ratio = phi


        self.solve_basic_properties()

        return None


    def solve_basic_properties(self):
        """
        """
        c = self.c
        f = self.far
        phi = self.equiv_ratio

        # Fuel flow:
        if not (c is None) and (f is None):
            c_ = c
            f_ = c_/self.mflow_e

        elif (c is None) and not (f is None):
            f_ = f
            c_ = f_*self.mflow_e
            
        elif not (c is None) and not (f is None):
            if np.abs(c/self.mflow_e-f)>1e-4:
                warnings.warn('Discrepancy between fuel flow ({})kg/s and FAR ({}). Fuel/Air ratio will be recalculated'.format(c, f))
                c_ = c
                f_ = c_/self.mflow_e

            else:
                c_ = c
                f_ = f
        
        if not (1e-3<f_<1e-1):
            warnings.warn('Unfeasible fuel/air ratio: {}'.format(f_))

        ## Store fuel flow and FAR:
        self._f = f_
        self._c = c_


        return


    def solve_combustion_chamber(self):
        """
        """

        return