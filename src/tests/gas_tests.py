"""
pyturb
Gas tests

M Rodriguez. 2020
"""

import numpy as np
import pyturb.utils.constants as cts
import pyturb.utils.units as units
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.perfect_ideal_gas import SemiperfectIdealGas
from pyturb.gas_models.thermo_properties import ThermoProperties

print('Perfect Gas ----------------------------------------------')

gas = 'Air'
tp = ThermoProperties(gas)

print(tp)

MO2 = 31.9988000 #g/mol
MN2 = 28.0134800 #g/mol
MAr = 39.9474514 #g/mol
MCO2 = 44.0095000 #g/mol
Ma = (78.084*MN2 + 20.9476*MO2 + 0.9365*MAr + 0.0319*MCO2 ) / 100

Ra = cts.Ru/Ma*1e3 # J/kg/K
print(Ra, Ma)

gamma_a = 1.4
cv = Ra/(gamma_a - 1)
cp = cv + Ra

pig = PerfectIdealGas(gas)
print(pig.gamma())
np.testing.assert_almost_equal(tp.Mg, Ma, decimal=4)
np.testing.assert_almost_equal(pig.Rg, Ra, decimal=4)
np.testing.assert_almost_equal(pig.cp(), cp, decimal=0)

spig = SemiperfectIdealGas(gas)

print(cv, spig.cv(cts.T_ref))

print(spig.cp(cts.T_ref), pig.cp())
print(spig.cv(cts.T_ref)*cts.T_ref, cv*cts.T_ref)
print(spig.gamma(cts.T_ref))

temp = 500 #K
h0 = cp*temp
h_pig = pig.cp()*temp
h_spig = spig.cp(temp)*temp

print(cp, pig.cp(), spig.cp(temp))
print(h0, h_pig, h_spig)


print('*******************************************')
h0_molar = spig.h0_molar(temp)

hf0Tref = tp.deltaHf_ref
h_molar_ref = spig.cp_molar(298.15)*298.15

H_molar_1 = hf0Tref + (h0_molar - h_molar_ref)
H_molar_2 = hf0Tref + (spig.cp_molar(temp)*temp - h_molar_ref)

print('H_molar_1=', H_molar_1)
print('H_molar_2=', H_molar_2)
print('h0_molar=', h0_molar)
print('h_molar_ref=', h_molar_ref)
print('hfoTref=', hf0Tref)


print('+++++++++++++++++++')
a1 = 1.009950160e+04
a2 = -1.968275610e+02
a3 =5.009155110e+0
a4 = -5.761013730e-03
a5 = 1.066859930e-05
a6 = -7.940297970e-09
a7 = 2.185231910e-12
b1 = -1.767967310e+02

temp = 500
T= temp
cp_ = (a1*T**-2 + a2 * T**-1 + a3 + a4 *T + a5 *T**2 + a6*T**3 + a7*T**4)*cts.Ru
hid = (cp_ * T)
h_ad = -a1*T**-2 + a2*np.log(T)/T + a3 + a4/2 * T + a5/3 * T**2 + a6/4 * T**3 + a7/5 * T**4 + b1/T
h_re = h_ad * cts.Ru * T

print(hid)
print(h_ad, spig.h0_dimensionless(temp))
print(h_re, spig.h0_molar(temp), pig.cp()*temp, spig.cp(temp)*temp)



