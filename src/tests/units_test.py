"""
pyturb
units.py tests

M Rodriguez
"""

import numpy as np
import pyturb.utils.units as units
import pyturb.utils.constants as cts

# Temperature
np.testing.assert_almost_equal(units.celsius_to_kelvin(15), 288.15)
np.testing.assert_almost_equal(units.kelvin_to_celsius(288.15), 15)
np.testing.assert_almost_equal(units.celsius_to_fahrenheit(20), 68)
np.testing.assert_almost_equal(units.fahrenheit_to_celsius(68), 20)

np.testing.assert_almost_equal(units.fahrenheit_to_kelvin(451), 505.928, decimal=-3)
np.testing.assert_almost_equal(units.kelvin_to_fahrenheit(373.15), 212)
print(units.celsius_to_fahrenheit(units.kelvin_to_celsius(373.15)))

np.testing.assert_almost_equal(units.celsius_to_rankine(15), 518.67)
np.testing.assert_almost_equal(units.rankine_to_celsius(1000), 282.40555555)
np.testing.assert_almost_equal(units.rankine_to_kelvin(1), 0.55555555)
np.testing.assert_almost_equal(units.kelvin_to_rankine(273.15), 536.67, decimal=-2)
print(units.fahrenheit_to_rankine(units.celsius_to_fahrenheit(0)))


# Pressure
np.testing.assert_almost_equal(units.Torr_to_Pa, 101325/760)
np.testing.assert_almost_equal(units.Pa_to_mmHg*cts.p_ref_SL, 760)


# 
