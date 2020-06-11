"""
"""

from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.isentropic_flow import IsentropicFlow

air = PerfectIdealGas('Air')
isent = IsentropicFlow(air)

#
T = 288.15 #K
p = 101325 #Pa
rho = p/air.Rg/T
v = 100 # m/s

a = isent.sound_speed(T)
M = isent.mach_number(v, T)
print('M=', M, 'a=', a)

print('Tt_T=', isent.stagnation_static_rel(M))
print('Tt=', isent.stag_temp_from_mach(M, T))
print('Tt=', isent.stag_temp_from_vel(v, T))
print('pt=', isent.stag_pressure_from_mach(M, p))
print('qi=', isent.impact_pressure_from_mach(M, p))
print('qi=', isent.dynamic_pressure(v, rho))
print('v=', isent.vel_from_stag_temp(isent.stag_temp_from_mach(M, T), T))
print('v=', isent.vel_from_mach(M, T))
print('rho=', rho, 'rhot=', isent.stag_density_from_mach(M, p/T/287))