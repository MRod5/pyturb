"""
pyturb
isa.py tests

M Rodriguez
"""

import numpy as np
import pyturb.units as units

np.testing.assert_almost_equal(units.celsius2kelvin(15), 288.15)
np.testing.assert_almost_equal(units.kelvin2celsius(288.15), 15)
np.testing.assert_almost_equal(units.celsius2fahrenheit(20), 68)
np.testing.assert_almost_equal(units.fahrenheit2celsius(68), 20)

np.testing.assert_almost_equal(units.fahrenheit2kelvin(212), 373.15)
np.testing.assert_almost_equal(units.kelvin2fahrenheit(373.15), units.celsius2fahrenheit(units.kelvin2celsius(373.15)))

np.testing.assert_almost_equal(units.celsius2rankine(15), 518.67)
np.testing.assert_almost_equal(units.rankine2celsius(1000), 282.40555555)
np.testing.assert_almost_equal(units.rankine2kelvin(1), 0.55555555)
np.testing.assert_almost_equal(units.kelvin2rankine(273.15), units.fahrenheit2rankine(units.celsius2fahrenheit(0)))

