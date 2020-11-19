"""
pyturb
Gas Mixture tests

M Rodriguez. 2020
"""

import pyturb.utils.constants as cts
from pyturb.gas_models import PerfectIdealGas
from pyturb.gas_models.gas_mixture import GasMixture
import numpy as np


def test_constants_are_equal():
    air = PerfectIdealGas('Air')
    gm = GasMixture(gas_model='Perfect')

    gm.add_gas('O2', 0.209476)
    gm.add_gas('N2', 0.78084)
    gm.add_gas('Ar', 0.009365)
    gm.add_gas('CO2', 0.000319)

    np.testing.assert_almost_equal(gm.Rg, air.Rg, decimal=6)
    np.testing.assert_almost_equal(gm.Mg, air.thermo_prop.Mg, decimal=6)
    np.testing.assert_almost_equal(gm.cp(cts.T_ref), air.cp(cts.T_ref), decimal=6)
    np.testing.assert_almost_equal(gm.cv(cts.T_ref), air.cv(cts.T_ref), decimal=6)
    np.testing.assert_almost_equal(gm.gamma(cts.T_ref), air.gamma(cts.T_ref), decimal=6)
