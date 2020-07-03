"""
"""
from pyturb.gas_models.gas_mixture import Gas
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
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
    

        # Mixture initialization:
        if not mixture is None:
            if not type(mixture) is dict:
                warnings.warn("mixture must be a dictionary with keys:=species and value:=moles. Instead received {}".format(mixture))

        return


    def add_gas(self, species, moles):
        """
        """
        return
