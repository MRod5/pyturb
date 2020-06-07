"""
Perfect Ideal Gas:
------------------

Perfect gas approach with Ideal Gas: cp, cv are constant and the Ideal Gas Law (pv=RgT) and
the Mayer Equation (cp-cv=Rg) are applicable:

    u = cv·T
    h = cp·T
    pv = Rg·T
    h = u + pv --> cp·T = cv·T + Rg·T --> Rg = cp - cv

MRodriguez. 2020

"""

import pyturb.utils.constants as cts
from pyturb.gas_models.gas import Gas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
from pyturb.gas_models.thermo_properties import ThermoProperties


class PerfectIdealGas(Gas):
    """
    PerfectIdealGas class:
    ----------------------

    Implements a Gas object:
        +gas_species: Species selected. May be a pure substance or any of the molecules and mixes
                      considered in "NASA Glenn Coefficients for Calculating Thermodynamic
                      Properties of Individual Species".
        + Ru: Ideal gas law universal constant  [J/mol/K]
        + Rg: Individual gas constant [J/kg/K]
        + cp: Constant heat capacity at constant pressure [J/kg/K]
        + cv: Constant heat capacity at constant volume [J/kg/K]
        + gamma: Heat capacity ratio [-]
    
    When initialized, SemiperfectIdealGas class is called to retrieve the gas species.
    The heat capacity at constant pressure is calculated as SemiperfectIdealGas.cp(T_ref)
    where T_ref is the reference temperatue for chemical processes.

    Heat capacity at constant volume and heat capacity ratio are obtained as an Ideal Gas:
    cv = cp - Rg, gamma=cp/cv.
    
    """

    def __init__(self, species=''):
        """
        Initializes a SemiperfectIdealGas object with the gas species selected and stores the
        gas species and the reference temperature.
        """

        # Check if the gas species is valid
        if species=='':
            raise ValueError("Gas species not defined.")
        if type(species) is not str:
            raise TypeError("Gas species must be a string")

        # SemiperfectIdealGas object:
        self.__from_semiperfect_gas = SemiperfectIdealGas(species)
        
        # Store the thermodynamic properties of the gas
        self.thermo_prop = self.__from_semiperfect_gas.thermo_prop

        # Reference temperature and gas species
        self.T_ref = cts.T_ref
        self.__gas_species = species
        return
        
    
    @property
    def gas_species(self):
        """
        Gets the Name of the gas species selected. May be a pure substance or any of the 
        molecules and mixes considered in "NASA Glenn Coefficients for Calculating Thermodynamic
        Properties of Individual Species".
        """
        return self.__gas_species
        

    @property
    def Ru(self):
        """
        Get the Ideal Gas Law constant Ru [J/mol/K]
        """
        Ru = cts.Ru
        return Ru
        
        
    @property
    def Rg(self):
        """
        Get the Individual Gas constant Rg =  Ru/Mg [J/kg/K]
        """
        Rg = self.Ru/self.thermo_prop.Mg*1e3
        return Rg
        
    
    def cp(self, temperature=None):
        """
        Heat capacity ratio at constant pressure [J/kg/K].

        As a perfect gas, cp is considered invariant with temperature. It is calculated as a
        semi-perfect gas (which means cp(T)) at T_ref temperature for any temperature.
        """
        cp_ = self.__from_semiperfect_gas.cp(self.T_ref)
        return cp_
        
    
    def cv(self, temperature=None):
        """
        Heat capacity ratio at constant volume [J/kg/K]

        As an ideal gas, cv is invariant with tempeature. cv is calculated with the 
        Mayer equation: cv = cp - Rg.
        """
        cv_ = self.cp(self.T_ref) - self.Rg
        return cv_
    
        
    def gamma(self, temperature=None):
        """
        Heat capacity ratio cp/cv [-].

        As an ideal gas, gamma is considered invariant with temperature. gamma is calculated
        as gamma = cp/cv.
        """
        gamma_ = self.cp(self.T_ref)/self.cv(self.T_ref)
        return gamma_
        
        