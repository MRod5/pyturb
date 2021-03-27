"""
pyturb

Tests for Semiperfect Gas model

M Rodriguez 2021
"""

import numpy as np
from pyturb.utils import constants as cts
from pyturb.gas_models import PerfectIdealGas
from pyturb.gas_models import SemiperfectIdealGas
from pyturb.gas_models import ThermoProperties

print('Perfect Gas ----------------------------------------------')

gas = 'Air'
tp = ThermoProperties(gas)
perfect_air = PerfectIdealGas(gas)
semiperfect_air = SemiperfectIdealGas(gas)


np.testing.assert_equal(perfect_air.Rg, semiperfect_air.Rg)

