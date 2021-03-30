"""
pyturb
Gas Mixture tests

M Rodriguez. 2020
"""

import pyturb.utils.constants as cts
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.perfect_ideal_gas import SemiperfectIdealGas
from pyturb.gas_models.gas_mixture import GasMixture
import numpy as np

air = PerfectIdealGas('Air')
gas_mixture_p_air = GasMixture(gas_model='Perfect')

gas_mixture_p_air.add_gas('O2', 0.209476)
gas_mixture_p_air.add_gas('N2', 0.78084)
gas_mixture_p_air.add_gas('Ar', 0.009365)
gas_mixture_p_air.add_gas('CO2', 0.000319)


print('--------------------------------------------------------------------------------------------------------------')
print(gas_mixture_p_air.mixture_gases)
print('--------------------------------------------------------------------------------------------------------------')

np.testing.assert_approx_equal(gas_mixture_p_air.Rg, air.Rg, significant=8)
np.testing.assert_approx_equal(gas_mixture_p_air.Mg, air.thermo_prop.Mg, significant=8)
np.testing.assert_approx_equal(gas_mixture_p_air.cp(cts.T_ref), air.cp(cts.T_ref), significant=8)
np.testing.assert_approx_equal(gas_mixture_p_air.cv(cts.T_ref), air.cv(cts.T_ref), significant=8)
np.testing.assert_approx_equal(gas_mixture_p_air.gamma(cts.T_ref), air.gamma(cts.T_ref), significant=8)

air = PerfectIdealGas('Air')
gas_mixture_sp_air = GasMixture(gas_model='Semi-perfect')

gas_mixture_sp_air.add_gas('O2', 0.209476)
gas_mixture_sp_air.add_gas('N2', 0.78084)
gas_mixture_sp_air.add_gas('Ar', 0.009365)
gas_mixture_sp_air.add_gas('CO2', 0.000319)


print('--------------------------------------------------------------------------------------------------------------')
print(gas_mixture_sp_air.mixture_gases)
print('--------------------------------------------------------------------------------------------------------------')

np.testing.assert_approx_equal(gas_mixture_sp_air.Rg, air.Rg, significant=8)
np.testing.assert_approx_equal(gas_mixture_sp_air.Mg, air.thermo_prop.Mg, significant=8)
np.testing.assert_approx_equal(gas_mixture_sp_air.cp(cts.T_ref), air.cp(cts.T_ref), significant=8)
np.testing.assert_approx_equal(gas_mixture_sp_air.cv(cts.T_ref), air.cv(cts.T_ref), significant=8)
np.testing.assert_approx_equal(gas_mixture_sp_air.gamma(cts.T_ref), air.gamma(cts.T_ref), significant=8)

print(gas_mixture_sp_air.h0(cts.T_ref))
print(air.h0(cts.T_ref))

semiperfect_air = SemiperfectIdealGas('Air')

print('---')

print(semiperfect_air.cp(280))
print(semiperfect_air.cp(1000))
print(semiperfect_air.cv(1500))
print(semiperfect_air.cv(2000))

np.testing.assert_approx_equal(semiperfect_air.cp(280), gas_mixture_sp_air.cp(280), significant=8)
np.testing.assert_approx_equal(semiperfect_air.cp(1000), gas_mixture_sp_air.cp(1000), significant=8)
np.testing.assert_approx_equal(semiperfect_air.cp(1500), gas_mixture_sp_air.cp(1500), significant=8)
np.testing.assert_approx_equal(semiperfect_air.cp(2000), gas_mixture_sp_air.cp(2000), significant=8)