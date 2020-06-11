"""
IsentropicFlow class:
---------------------

Isentropic flow relations


MRodriguez 2020

"""

from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
import numpy as np


class IsentropicFlow(object):
    """
    Isentropic flow relations.

    """

    def __init__(self, fluid):
        if not(isinstance(fluid, PerfectIdealGas) or isinstance(fluid, SemiperfectIdealGas)):
            raise TypeError("Object must be PerfectIdealGas, SemiperfectIdealGas. Instead received {}".format(fluid))
        
        for attr in ['cp', 'gamma']:
            if not hasattr(fluid, attr):
                raise AttributeError("Attribute {} not found in fluid object".format(attr))

        self.fluid = fluid
        return