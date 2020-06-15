from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.power_plant.intake import Intake
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

print('h0={0:6.2f}m, v0={1:5.3f}m/s, p0={2:8.1f}Pa, T0={3:5.1f}K, rho0={4:6.4f}kg/m**3'.format(h0, v0, p0, T0, intake.rho_e))
print('M1={0:5.3f}, p1t={1:8.1f}Pa, T1t={2:5.1f}K, rho1t={3:6.4f}kg/m**3'.format(intake.mach_e, intake.p_et, intake.T_et, intake.rho_et))
print('M2={0:5.3f}, p2t={1:8.1f}Pa, T2t={2:5.1f}K, rho2t={3:6.4f}kg/m**3'.format(intake.mach_s, intake.p_st, intake.T_st, intake.rho_st))
print('h1t={0:8.4f}W/kg/s, h1={1:8.4f}W/kg/s, ekin_1={2:8.4f}W/kg/s'.format(intake.h_et*1e-6, intake.h_e*1e-6, intake.ekin_e*1e-6))
print('h2t={0:8.4f}W/kg/s, h2={1:8.4f}W/kg/s, ekin_2={2:8.4f}W/kg/s'.format(intake.h_st*1e-6, intake.h_s*1e-6, intake.ekin_s*1e-6))
print('G1={0:5.1f}kg/s, G2={1:5.2f}kg/s, diffG21={2:5.2f}kg/s'.format(intake.mflow_e, intake.mflow_s, intake.delta_massflow))

