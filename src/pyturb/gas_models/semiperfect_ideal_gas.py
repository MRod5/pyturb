"""
"""

import numpy as np
import pyturb.utils.constants as cts
from pyturb.gas_models.gas import Gas
from pyturb.gas_models.thermo_properties import ThermoProperties


class SemiperfectIdealGas(Gas):
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
        tmin = np.min(self.thermo_prop.temp_range)
        tmax = np.max(self.thermo_prop.temp_range)
        if tmin>temperature or tmax<temperature:
            raise ValueError("Gas temperature ({0}K) out of implemented limits [{1},{2}]K".format(temperature, tmin, tmax))

        for (temp_range, coeffs) in zip(self.thermo_prop.temp_range, self.thermo_prop.coefficients):
            if temp_range[0]<=temperature<temp_range[1]:
                temp_poly = np.array([temperature**(-2), temperature**(-1), 1, temperature, temperature**(2), temperature**(3), temperature**(4)])
                cp_ = np.dot(coeffs, temp_poly)
                
                cp_ = cp_ * self.Rg
        return cp_
        
    
    def cv(self, temperature):
        cv_ = self.cp(temperature) - self.Rg
        return cv_
    
        
    def gamma(self, temperature):
        gamma_ = self.cp(temperature)/self.cv(temperature)
        return gamma_
        
        