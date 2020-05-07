"""
pyturb
isa.py tests

M Rodriguez
"""

import numpy as np
import pyturb.isa as isa

np.testing.assert_almost_equal(isa.temp_isa(0), [15+273.15])
np.testing.assert_almost_equal(isa.pres_isa(0),[101325])
np.testing.assert_almost_equal(isa.temp_isa(11000), [216.65])
np.testing.assert_almost_equal(isa.pres_isa(11000), 22700.0, decimal=-2)

