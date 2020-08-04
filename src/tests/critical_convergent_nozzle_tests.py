"""
"""

from pyturb.gas_models.isentropic_flow import IsentropicFlow
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.power_plant.nozzle import Nozzle
import pyturb.gas_models.isa as isa
import numpy as np

air = PerfectIdealGas('Air')
isent_flow = IsentropicFlow(air)

h0 = 0.0
p0 = isa.pressure_isa(h0)
T0 = isa.temperature_isa(h0)
phi_1 = 1 # m
Ae = np.pi * phi_1**2/4
print(Ae)
v0 = 300
G0 = p0/T0/air.Rg * Ae * v0
M0 = isent_flow.mach_number(v0, T0)
p0t = isent_flow.stag_pressure_from_mach(M0, p0, T0)
T0t = isent_flow.stag_temp_from_mach(M0, T0)

nozzle = Nozzle(air)
nozzle.initialize_nozzle(G0, p0t, T0t, Ae=Ae)
nozzle.solve_critical_convergent_nozzle(adiab_efficiency=0.9)


print('--- CRITICAL CONVERGENT NOZZLE: calc ps, As ----------------------------------------------------------------------')
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

nozzle.solve_critical_convergent_nozzle(ps=75000)
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

nozzle.solve_critical_convergent_nozzle(As=0.9)
print('--- CRITICAL CONVERGENT NOZZLE: calc As, adiab_eff ----------------------------------------------------------------------')
print('M7={0:5.3f}, p7t={1:8.1f}Pa, T7t={2:5.1f}K, rho7t={3:6.4f}kg/m**3'.format(nozzle.mach_e, nozzle.p_et, nozzle.T_et, nozzle.rho_et))
print('v7={0:5.2f}, p7={1:8.1f}Pa, T7={2:5.1f}K, rho7={3:6.4f}kg/m**3'.format(nozzle.vel_e, nozzle.p_e, nozzle.T_e, nozzle.rho_e))
print('M8={0:5.3f}, p8t={1:8.1f}Pa, T8t={2:5.1f}K, rho8t={3:6.4f}kg/m**3'.format(nozzle.mach_s, nozzle.p_st, nozzle.T_st, nozzle.rho_st))
print('v8={0:5.2f}, p8={1:8.1f}Pa, T8={2:5.1f}K, rho8={3:6.4f}kg/m**3'.format(nozzle.vel_s, nozzle.p_s, nozzle.T_s, nozzle.rho_s))
print('h7t={0:8.4f}kW/kg/s, h7={1:8.4f}kW/kg/s, ekin_7={2:8.4f}kW/kg/s'.format(nozzle.h_et*1e-3, nozzle.h_e*1e-3, nozzle.ekin_e*1e-3))
print('h8t={0:8.4f}kW/kg/s, h8={1:8.4f}kW/kg/s, ekin_8={2:8.4f}kW/kg/s'.format(nozzle.h_st*1e-3, nozzle.h_s*1e-3, nozzle.ekin_s*1e-3))
print('G7={0:5.1f}kg/s, G8={1:5.2f}kg/s, diffG87={2:5.2f}kg/s'.format(nozzle.mflow_e, nozzle.mflow_s, nozzle.delta_massflow))
print('A9={0:8.4f}m**2, adiab_eff={1}'.format(nozzle.A_s, nozzle.adiab_efficiency))
print('p_star={0:8.1f}Pa, T_star={1:5.1f}, A_star={2:8.4f}'.format(nozzle.p_s_star, nozzle.T_s_star, nozzle.A_star))

nozzle.solve_critical_convergent_nozzle(As=np.pi, ps=np.pi)
print('--- CRITICAL CONVERGENT NOZZLE: calc As, adiab_eff ----------------------------------------------------------------------')
print('M7={0:5.3f}, p7t={1:8.1f}Pa, T7t={2:5.1f}K, rho7t={3:6.4f}kg/m**3'.format(nozzle.mach_e, nozzle.p_et, nozzle.T_et, nozzle.rho_et))
print('v7={0:5.2f}, p7={1:8.1f}Pa, T7={2:5.1f}K, rho7={3:6.4f}kg/m**3'.format(nozzle.vel_e, nozzle.p_e, nozzle.T_e, nozzle.rho_e))
print('M8={0:5.3f}, p8t={1:8.1f}Pa, T8t={2:5.1f}K, rho8t={3:6.4f}kg/m**3'.format(nozzle.mach_s, nozzle.p_st, nozzle.T_st, nozzle.rho_st))
print('v8={0:5.2f}, p8={1:8.1f}Pa, T8={2:5.1f}K, rho8={3:6.4f}kg/m**3'.format(nozzle.vel_s, nozzle.p_s, nozzle.T_s, nozzle.rho_s))
print('h7t={0:8.4f}kW/kg/s, h7={1:8.4f}kW/kg/s, ekin_7={2:8.4f}kW/kg/s'.format(nozzle.h_et*1e-3, nozzle.h_e*1e-3, nozzle.ekin_e*1e-3))
print('h8t={0:8.4f}kW/kg/s, h8={1:8.4f}kW/kg/s, ekin_8={2:8.4f}kW/kg/s'.format(nozzle.h_st*1e-3, nozzle.h_s*1e-3, nozzle.ekin_s*1e-3))
print('G7={0:5.1f}kg/s, G8={1:5.2f}kg/s, diffG87={2:5.2f}kg/s'.format(nozzle.mflow_e, nozzle.mflow_s, nozzle.delta_massflow))
print('A9={0:8.4f}m**2, adiab_eff={1}'.format(nozzle.A_s, nozzle.adiab_efficiency))
print('p_star={0:8.1f}Pa, T_star={1:5.1f}, A_star={2:8.4f}'.format(nozzle.p_s_star, nozzle.T_s_star, nozzle.A_star))

