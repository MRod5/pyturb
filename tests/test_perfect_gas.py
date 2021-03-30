"""
pyturb

Tests for Perfect Gas model

M Rodriguez 2021
"""

import numpy as np
from pyturb.utils import constants as cts
from pyturb.gas_models import PerfectIdealGas
from pyturb.gas_models import ThermoProperties

print('Perfect Gas ----------------------------------------------')

gas = 'Air'
tp = ThermoProperties(gas)
pig = PerfectIdealGas(gas)

print(tp)

# From NASA Earth fact sheet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
MO2 = 31.9988000 #g/mol
MN2 = 28.0134800 #g/mol
MAr = 39.9474514 #g/mol
MCO2 = 44.0095000 #g/mol
MNe = 20.1797 #g/mol
Ma = (78.084*MN2 + 20.9476*MO2 + 0.9365*MAr + 0.0319*MCO2 ) / 100

print((78.084 + 20.9476 + 0.9365 + 0.0319 ))
Ra = cts.Ru/Ma*1e3 # J/kg/K

print(Ra, Ma)

gamma_a = pig.gamma()
cv = Ra/(gamma_a - 1)
cp = cv + Ra

print(pig.gamma(), pig.Ru)
np.testing.assert_approx_equal(pig.Ru, 8.31446261815324, significant=12)
np.testing.assert_approx_equal(tp.Mg, Ma, significant=7)
np.testing.assert_approx_equal(pig.Rg, Ra, significant=7)
np.testing.assert_approx_equal(pig.cp(), cp, significant=7)

print(cv, pig.cv())
print(cp, pig.cp())

print('*******************************************')


