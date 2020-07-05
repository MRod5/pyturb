"""
"""
import pyturb.utils.constants as cts
from pyturb.gas_models.gas import Gas
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
import numpy as np
import pandas as pd
import warnings


class GasMixture(object):
    """
    """

    def __init__(self, gas_model = "perfect", mixture=None):
        """
        """

        # Gas model selector:
        if gas_model.lower() == "perfect":
            self.gas_model = PerfectIdealGas
        elif gas_model.lower() == "semi-perfect" or gas_model.lower() == "semiperfect":
            self.gas_model = SemiperfectIdealGas
        else:
            self.gas_model = None
            raise ValueError("gas_model may be 'perfect' or 'semi-perfect', instead received {}".format(gas_model))
        
        columns = ['pure_substance_props', 'Ng', 'Mg', 'Rg']
        self.mixture_gases = pd.DataFrame(columns=columns, index=['gas0'])

        # Mixture initialization:
        if not mixture is None:
            if not type(mixture) is dict:
                warnings.warn("mixture must be a dictionary with keys:=species and value:=moles. Instead received {}".format(mixture))
            else:
                #TODO: call add_gas for each pure substance in the mixture
                print("")


        return


    def add_gas(self, species, moles):
        """
        """

        pure_substance = self.gas_model(species)

        # TODO: Evaluate using pandas instead...
        self.mixture_gases.append({'substance_props': pure_substance,'Ng': moles, 'Mg': pure_substance.thermo_prop.Mg, 'Rg': pure_substance.Rg}, ignore_index=True)
        print(self.mixture_gases)
        #self._get_mixture_properties()

        return


    def _get_mixture_properties(self):
        """
        """
        self._Nmix = 0
        self._Mg = 0
        
        for gas in self.pure_substance:
            pure_substance = self.pure_substance[gas]
            self._Nmix += pure_substance['moles']
            # TODO: dataframe will be more efficient here, one for-loop would be needed
            
        
        for gas in self.pure_substance:
            pure_substance = self.pure_substance[gas]
            molar_frac[gas] = pure_substance['moles']/self.Nmix
            self._Mg += pure_substance['substance_props'].thermo_prop.Mg * molar_frac[gas]


        return


    @property
    def Nmix(self):
        """
        """
        return self._Nmix
    

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
        Get the Mixture Gas constant Rg =  Ru/Mg [J/kg/K]
        """
        
        Rg = self.Ru/self.Mg*1e3
        return Rg


    @property
    def Mg(self):
        """
        Get the Mixture molecular mass [g/mol]
        """
        return self._Mg