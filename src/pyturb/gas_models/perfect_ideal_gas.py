"""
"""

import pyturb.utils.constants as cts
from pyturb.gas_models.gas import Gas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
from pyturb.gas_models.thermo_properties import ThermoProperties


class PerfectIdealGas(Gas):
    """
    """

    def __init__(self, species=''):
        """
        """

        if species=='':
            raise ValueError("Gas species not defined.")
        if type(species) is not str:
            raise TypeError("Gas species must be a string")


        # self.thermo_prop = ThermoProperties(species)
        # if self.thermo_prop.Mg == None:
        #     raise ValueError("Gas unavailable: {}".format(species))
        # else:
        #     self.__gas_species = species

        self.__from_semiperfect_gas = SemiperfectIdealGas(species)
        self.thermo_prop = self.__from_semiperfect_gas.thermo_prop
        self.T_ref = cts.T_ref
        self.__gas_species = species
        return
        
    
    @property
    def gas_species(self):
        return self.__gas_species
        

    @property
    def Ru(self):
        Ru = cts.Ru
        return Ru
        
        
    @property
    def Rg(self):
        Rg = self.Ru/self.thermo_prop.Mg*1e3
        return Rg
        
    
    def cp(self, temperature=None):
        cp_ = self.__from_semiperfect_gas.cp(self.T_ref)
        return cp_
        
    
    def cv(self, temperature=None):
        cv_ = self.cp(self.T_ref) - self.Rg
        return cv_
    
        
    def gamma(self, temperature=None):
        gamma_ = self.cp(self.T_ref)/self.cv(self.T_ref)
        return gamma_
        
        