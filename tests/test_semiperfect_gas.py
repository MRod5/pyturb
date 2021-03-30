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
perfect_air = PerfectIdealGas(gas)
semiperfect_air = SemiperfectIdealGas(gas)

print(semiperfect_air.thermo_prop)

np.testing.assert_equal(perfect_air.gas_species, semiperfect_air.gas_species)
np.testing.assert_equal(perfect_air.Rg, semiperfect_air.Rg)
np.testing.assert_equal(perfect_air.Mg, semiperfect_air.Mg)
np.testing.assert_equal(perfect_air.cp(), semiperfect_air.cp(cts.T_ref))

np.testing.assert_approx_equal(semiperfect_air.cp_molar(cts.T_ref)/semiperfect_air.Mg*1e3, perfect_air.cp(),significant=10)
np.testing.assert_approx_equal(semiperfect_air.cv_molar(cts.T_ref)/semiperfect_air.Mg*1e3, perfect_air.cv(),significant=10)

np.testing.assert_equal(perfect_air.gamma(), semiperfect_air.gamma(cts.T_ref))

print(semiperfect_air.cp(1000))
print(semiperfect_air.cp(2000))
print(semiperfect_air.cv(1000))
print(semiperfect_air.cv(2000))

print(semiperfect_air.cp(1000)-semiperfect_air.cv(1000))
print(semiperfect_air.cp(2000)-semiperfect_air.cv(2000))

np.testing.assert_approx_equal(semiperfect_air.cp(1000)-semiperfect_air.cv(1000), semiperfect_air.Rg,significant=10)
np.testing.assert_approx_equal(semiperfect_air.cp(2000)-semiperfect_air.cv(2000), semiperfect_air.Rg,significant=10)


temperature = np.array([200, cts.T_ref, cts.T_ref_SL, 1000, 1500, 2000, 2500])
print(semiperfect_air.cp(temperature))