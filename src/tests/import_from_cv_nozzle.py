"""
pyturb
import from other CV (Nozzle) tests

M Rodriguez. 2020
"""

import sys
from sys import path
from os.path import dirname as dir
sys.path.append(dir(sys.path[0]))

from pyturb.gas_models.isentropic_flow import IsentropicFlow
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.power_plant.nozzle import Nozzle
from pyturb.power_plant.nozzle import Intake
import pyturb.gas_models.isa as isa
import numpy as np

air = PerfectIdealGas('Air')
isent_flow = IsentropicFlow(air)

h0 = 5000.0 # m
p0 = isa.pressure_isa(h0)
T0 = isa.temperature_isa(h0)
phi_1 = 0.8 # m
Ae = np.pi * phi_1**2/4
v0 = 300 # m/s
G0 = p0/T0/air.Rg * Ae * v0
M0 = isent_flow.mach_number(v0, T0)
p0t = isent_flow.stag_pressure_from_mach(M0, p0, T0)
T0t = isent_flow.stag_temp_from_mach(M0, T0)


intake = Intake(air)
intake.initialize_intake(p0, T0, v0, Ae, As=0.75*Ae)


nozzle = Nozzle(air)
nozzle.initialize_nozzle(G0, p0t, T0t, Ae=Ae)
nozzle.solve_critical_convergent_nozzle(ps=50000)

Me = nozzle.mach_e
p7t = nozzle.p_et
T7t = nozzle.T_et
G7 = nozzle.mflow_s
A9 = nozzle.A_s
p9 = nozzle.p_s
T9 = nozzle.T_s

print('--- CRITICAL CONVERGENT NOZZLE: calc ps, adiab_eff ----------------------------------------------------------------------')
print('M7={0:5.3f}, p7t={1:8.1f}Pa, T7t={2:5.1f}K, rho7t={3:6.4f}kg/m**3'.format(nozzle.mach_e, nozzle.p_et, nozzle.T_et, nozzle.rho_et))
print('v7={0:5.2f}, p7={1:8.1f}Pa, T7={2:5.1f}K, rho7={3:6.4f}kg/m**3'.format(nozzle.vel_e, nozzle.p_e, nozzle.T_e, nozzle.rho_e))
print('M8={0:5.3f}, p8t={1:8.1f}Pa, T8t={2:5.1f}K, rho8t={3:6.4f}kg/m**3'.format(nozzle.mach_s, nozzle.p_st, nozzle.T_st, nozzle.rho_st))
print('v8={0:5.2f}, p8={1:8.1f}Pa, T8={2:5.1f}K, rho8={3:6.4f}kg/m**3'.format(nozzle.vel_s, nozzle.p_s, nozzle.T_s, nozzle.rho_s))
print('h7t={0:8.4f}kW/kg/s, h7={1:8.4f}kW/kg/s, ekin_7={2:8.4f}kW/kg/s'.format(nozzle.h_et*1e-3, nozzle.h_e*1e-3, nozzle.ekin_e*1e-3))
print('h8t={0:8.4f}kW/kg/s, h8={1:8.4f}kW/kg/s, ekin_8={2:8.4f}kW/kg/s'.format(nozzle.h_st*1e-3, nozzle.h_s*1e-3, nozzle.ekin_s*1e-3))
print('G7={0:5.1f}kg/s, G8={1:5.2f}kg/s, diffG87={2:5.2f}kg/s'.format(nozzle.mflow_e, nozzle.mflow_s, nozzle.delta_massflow))
print('A9={0:8.4f}m**2, adiab_eff={1}'.format(nozzle.A_s, nozzle.adiab_efficiency))
print('p_star={0:8.1f}Pa, T_star={1:5.1f}, A_star={2:8.4f}'.format(nozzle.p_s_star, nozzle.T_s_star, nozzle.A_star))
print('')

nozzle2 = Nozzle(air)
nozzle2.initialize_from_cv(intake)
nozzle2.solve_critical_convergent_nozzle(ps=50000)

np.testing.assert_almost_equal(nozzle2.mach_e, Me)
np.testing.assert_almost_equal(nozzle2.p_et, p7t)
np.testing.assert_almost_equal(nozzle2.T_et, T7t)
np.testing.assert_almost_equal(nozzle2.A_s, A9)
np.testing.assert_almost_equal(nozzle2.p_s, p9)
np.testing.assert_almost_equal(nozzle2.T_s, T9)