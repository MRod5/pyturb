"""
pyturb
Gas Mixture tests

M Rodriguez. 2020
"""

import pyturb.utils.constants as cts
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.gas_mixture import GasMixture
import numpy as np

air = PerfectIdealGas('Air')
gm = GasMixture(gas_model='Perfect')

gm.add_gas('O2', 0.209476)
gm.add_gas('N2', 0.78084)
gm.add_gas('Ar', 0.009365)
gm.add_gas('CO2', 0.000319)


print('--------------------------------------------------------------------------------------------------------------')
print(gm.mixture_gases)
print('--------------------------------------------------------------------------------------------------------------')

np.testing.assert_almost_equal(gm.Rg, air.Rg, decimal=6)
np.testing.assert_almost_equal(gm.Mg, air.thermo_prop.Mg, decimal=6)
np.testing.assert_almost_equal(gm.cp(cts.T_ref), air.cp(cts.T_ref), decimal=6)
np.testing.assert_almost_equal(gm.cv(cts.T_ref), air.cv(cts.T_ref), decimal=6)
np.testing.assert_almost_equal(gm.gamma(cts.T_ref), air.gamma(cts.T_ref), decimal=6)

air = PerfectIdealGas('Air')
gm2 = GasMixture(gas_model='Semi-perfect')

gm2.add_gas('O2', 0.209476)
gm2.add_gas('N2', 0.78084)
gm2.add_gas('Ar', 0.009365)
gm2.add_gas('CO2', 0.000319)


print('--------------------------------------------------------------------------------------------------------------')
print(gm2.mixture_gases)
print('--------------------------------------------------------------------------------------------------------------')

np.testing.assert_almost_equal(gm2.Rg, air.Rg, decimal=6)
np.testing.assert_almost_equal(gm2.Mg, air.thermo_prop.Mg, decimal=6)
np.testing.assert_almost_equal(gm2.cp(cts.T_ref), air.cp(cts.T_ref), decimal=6)
np.testing.assert_almost_equal(gm2.cv(cts.T_ref), air.cv(cts.T_ref), decimal=6)
np.testing.assert_almost_equal(gm2.gamma(cts.T_ref), air.gamma(cts.T_ref), decimal=6)

print(gm2.h0(cts.T_ref))
print(air.h0(cts.T_ref))
