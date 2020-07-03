"""
"""
import pyturb.utils.constants as cts
from pyturb.gas_models.gas_mixture import Gas
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
import numpy as np
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
        

        self.pure_substance = {}

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

        self.pure_substance[pure_substance.gas_species] = {'subtance_props': pure_substance,'moles': moles}

        self._get_mixture_properties()

        return


    def _get_mixture_properties(self):
        """
        """
        self._Nmix = np.sum(self.pure_subtance.value())
        return

    
    @property
    def Ru(self):
        """
        Get the Ideal Gas Law constant Ru [J/mol/K]
        """
        
        Ru = cts.Ru
        return Ru

    
    @property
    def Nmix(self):
        """
        """
        return self._Nmix