"""
pyturb
isa.py tests

M Rodriguez
"""

from isentropic_gas import IsentropicGas
import numpy as np

isentropic = IsentropicGas(selected_cp_air_model='ideal', selected_gamma_air_model='aaa')
isentropic = IsentropicGas(selected_cp_air_model='ideal', selected_gamma_air_model='ideal')
print('...........')
isentropic = IsentropicGas(selected_cp_air_model='nasa', selected_gamma_air_model='naca')

np.testing.assert_almost_equal(isentropic.stagnation_static_temperature_relation(0.5, 273), 1.05)
np.testing.assert_almost_equal(isentropic.stagnation_static_pressure_relation(0.5, 273), 1.186, decimal=3)

np.testing.assert_almost_equal(isentropic.kinetic_energy_from_enthalpy(15000, 10000), 5000, decimal=1)
np.testing.assert_almost_equal(isentropic.velocity_from_enthalpy(15000, 10000), 100, decimal=1)

np.testing.assert_almost_equal(isentropic.stagnation_enthalpy(100, 273), 279600, decimal=-3)
np.testing.assert_almost_equal(isentropic.stagnation_pressure_from_mach(0.5, 101325, 273), 101325*1.186, decimal=-3)
np.testing.assert_almost_equal(isentropic.mach_number(165.75, 273), 0.5, decimal=3)
np.testing.assert_almost_equal(isentropic.stagnation_temperature_from_mach(0.5, 273), 273*1.05, decimal=4)
np.testing.assert_almost_equal(isentropic.stagnation_temperature_from_vel(165.75, 273), 273*1.05, decimal=2)
