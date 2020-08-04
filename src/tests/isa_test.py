"""
pyturb
isa.py tests

M Rodriguez 2020
"""

import numpy as np
import pyturb.gas_models.isa  as isa

isa.temperature_isa([0, 1000, 2000])


h = np.linspace(0,85000,100)
T = isa.temperature_isa(h)
p = isa.pressure_isa(h)
d = isa.density_state_eq(h)


h = isa.height_from_temperature_isa(288.15)
h2 = isa.height_from_pressure_isa(85000)
hh = isa.height_from_pressure_isa([50000, 22000, 1500])
print(h, h2)
isa.height_from_temperature_isa(270)
isa.height_from_temperature_isa([220, 250])


np.testing.assert_almost_equal(isa.temperature_isa(0), [15+273.15])
np.testing.assert_almost_equal(isa.pressure_isa(0),[101325])
np.testing.assert_almost_equal(isa.temperature_isa(11000), [216.65])
np.testing.assert_almost_equal(isa.pressure_isa(11000), 22700.0, decimal=-2)


T75 = isa.temperature_isa(7500)
T90 = isa.temperature_isa(9000, 10)
p75 = isa.pressure_isa(7500)
p90 = isa.pressure_isa(9000,10)

print(T75, T90)
print(p75, p90)

hp75 = isa.height_from_pressure_isa(p75)
hp90 = isa.height_from_pressure_isa(p90, 10)
print(hp75, hp90)

np.testing.assert_almost_equal(hp75,7500)
np.testing.assert_almost_equal(hp90,9000)