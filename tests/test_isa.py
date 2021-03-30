"""
pyturb

Tests for isa.py

M Rodriguez 2021
"""

import numpy as np
from pyturb.gas_models import isa


isa.temperature_isa([0, 1000, 2000])

h = np.linspace(0,11000,21)
T = isa.temperature_isa(h)
p = isa.pressure_isa(h)
d = isa.density_isa(h)

np.testing.assert_almost_equal(d[0], 1.2250, decimal=4)
np.testing.assert_almost_equal(d[-1], 0.3639, decimal=4)
np.testing.assert_almost_equal(p[-1], 22632)

h0 = isa.height_from_temperature_isa(288.15)
h11 = isa.height_from_pressure_isa(22632)
h32 = isa.height_from_pressure_isa(868.02)
hh = isa.height_from_pressure_isa([110.91, 66.939, 3.9564])

np.testing.assert_almost_equal(h0, 0.0)
np.testing.assert_almost_equal(h11, 11000.0)
np.testing.assert_array_almost_equal(hh, [47000, 51000, 71000], decimal=5)
np.testing.assert_almost_equal(isa.pressure_isa(0),[101325])

np.testing.assert_almost_equal(isa.temperature_isa(0), [15+273.15])
np.testing.assert_almost_equal(isa.temperature_isa(11000), [216.65])
np.testing.assert_array_almost_equal(h, isa.height_from_pressure_isa(p))


T75 = isa.temperature_isa(7500)
T90 = isa.temperature_isa(9000, 9.75)
p75 = isa.pressure_isa(7500)
p90 = isa.pressure_isa(9000,9.75)

print(T75, T90)
print(p75, p90)

hp75 = isa.height_from_pressure_isa(p75)
hp90 = isa.height_from_pressure_isa(p90, 9.75)
print(hp75, hp90)

np.testing.assert_almost_equal(hp75,7500)
np.testing.assert_almost_equal(hp90,9000)

T = np.arange(273.15, 240, -5)
h = isa.height_from_temperature_isa(T)
print(h)