"""
"""

import pyturb.utils.constants as cts
from pyturb.gas_models.gas import Gas
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


        self.thermo_prop = ThermoProperties(species)
        if self.thermo_prop.Mg == None:
            raise ValueError("Gas unavailable: {}".format(species))
        else:
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
        
    
    def cp(self, temperature):
        return 1004
        
    
    def cv(self, temperature):
        cv_ = self.cp(temperature) - self.Rg
        return cv_
    
        
    def gamma(self, temperature):
        gamma_ = self.cp(temperature)/self.cv(temperature)
        return gamma_
        
        