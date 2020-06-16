from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.power_plant.intake import Intake
from pyturb.power_plant.nozzle import Nozzle
import pyturb.gas_models.isa as isa
import numpy as np


h0 = 0.0
p0 = isa.pressure_isa(h0)
T0 = isa.temperature_isa(h0)
phi_1 = 1 # m
Ae = np.pi * phi_1**2/4
As = 1.1*Ae
v0 = 300

air = PerfectIdealGas('Air')

intake = Intake(air)
intake.initialize_intake(p0, T0, v0, Ae, As=As)
print('--- INTAKE ------------------------------------------------------------------------------')
print('h0={0:6.2f}m, v0={1:5.3f}m/s, p0={2:8.1f}Pa, T0={3:5.1f}K, rho0={4:6.4f}kg/m**3'.format(h0, v0, p0, T0, intake.rho_e))
print('M1={0:5.3f}, p1t={1:8.1f}Pa, T1t={2:5.1f}K, rho1t={3:6.4f}kg/m**3'.format(intake.mach_e, intake.p_et, intake.T_et, intake.rho_et))
print('M2={0:5.3f}, p2t={1:8.1f}Pa, T2t={2:5.1f}K, rho2t={3:6.4f}kg/m**3'.format(intake.mach_s, intake.p_st, intake.T_st, intake.rho_st))
print('h1t={0:8.4f}W/kg/s, h1={1:8.4f}W/kg/s, ekin_1={2:8.4f}W/kg/s'.format(intake.h_et*1e-6, intake.h_e*1e-6, intake.ekin_e*1e-6))
print('h2t={0:8.4f}W/kg/s, h2={1:8.4f}W/kg/s, ekin_2={2:8.4f}W/kg/s'.format(intake.h_st*1e-6, intake.h_s*1e-6, intake.ekin_s*1e-6))
print('G1={0:5.1f}kg/s, G2={1:5.2f}kg/s, diffG21={2:5.2f}kg/s'.format(intake.mflow_e, intake.mflow_s, intake.delta_massflow))


G = intake.mflow_s
p7t = intake.p_st
T7t = intake.T_st

nozzle = Nozzle(air)
nozzle.initialize_nozzle(G, p7t, T7t)
nozzle.solve_adapted_nozzle(ps=p0, As=None)
print('--- NOZZLE ------------------------------------------------------------------------------')
print('M7={0:5.3f}, p7t={1:8.1f}Pa, T7t={2:5.1f}K, rho7t={3:6.4f}kg/m**3'.format(nozzle.mach_e, nozzle.p_et, nozzle.T_et, nozzle.rho_et))
print('M8={0:5.3f}, p8t={1:8.1f}Pa, T8t={2:5.1f}K, rho8t={3:6.4f}kg/m**3'.format(nozzle.mach_s, nozzle.p_st, nozzle.T_st, nozzle.rho_st))
print('v8={0:5.2f}, p8={1:8.1f}Pa, T8={2:5.1f}K, rho8={3:6.4f}kg/m**3'.format(nozzle.vel_s, nozzle.p_s, nozzle.T_s, nozzle.rho_s))
print('h7t={0:8.4f}W/kg/s, h7={1:8.4f}W/kg/s, ekin_7={2:8.4f}W/kg/s'.format(nozzle.h_et*1e-6, nozzle.h_e*1e-6, nozzle.ekin_e*1e-6))
print('h8t={0:8.4f}W/kg/s, h8={1:8.4f}W/kg/s, ekin_8={2:8.4f}W/kg/s'.format(nozzle.h_st*1e-6, nozzle.h_s*1e-6, nozzle.ekin_s*1e-6))
print('G7={0:5.1f}kg/s, G8={1:5.2f}kg/s, diffG87={2:5.2f}kg/s'.format(nozzle.mflow_e, nozzle.mflow_s, nozzle.delta_massflow))
print('A9={0:8.4f}m**2'.format(nozzle.A_s))
print('p_star={0:8.1f}Pa, T_star={1:5.1f}, A_star={2:8.4f}'.format(nozzle.p_s_star, nozzle.T_s_star, nozzle.A_star))

