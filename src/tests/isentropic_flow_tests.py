"""
pyturb
isentropic_flow tests

M Rodriguez. 2020
"""

from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas, SemiperfectIdealGas
from pyturb.gas_models.isentropic_flow import IsentropicFlow
import pyturb.gas_models.isa as isa

air = PerfectIdealGas('Air')
isent = IsentropicFlow(air)


# Perfect gas
print('----------------------------PerfectGas----------------------------')
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
print('rho=', rho, 'rhot=', isent.stag_density_from_mach(M, rho))
print('T=', isent.stat_temp_from_mach(M, isent.stag_temp_from_mach(M, T)))
print('p=', isent.stat_pressure_from_mach(M, isent.stag_pressure_from_mach(M, p)))

## Semiperfect gas
print('--------------------------SemiperfectGas--------------------------')
air = SemiperfectIdealGas('Air')
isent = IsentropicFlow(air)


##
h = 15000
T = isa.temperature_isa(h)
p = isa.pressure_isa(h)
rho = p/air.Rg/T
v = 500 # m/s

a = isent.sound_speed(T)
M = isent.mach_number(v, T)
print('M=', M, 'a=', a, 'T=', T, 'p=', p)

Tt_T = isent.stagnation_static_rel(M, T)
Tt = isent.stag_temp_from_mach(M, T)
pt = isent.stag_pressure_from_mach(M, p, T)
print('Tt_T=', Tt_T)
print('Tt=', Tt)
print('Tt=', isent.stag_temp_from_vel(v, T))
print('pt=', pt)
print('qi=', isent.impact_pressure_from_mach(M, p, T))
print('qi=', isent.dynamic_pressure(v, rho))
print('v=', isent.vel_from_stag_temp(isent.stag_temp_from_mach(M, T), T))
print('v=', isent.vel_from_mach(M, T))
print('rho=', rho, 'rhot=', isent.stag_density_from_mach(M, rho, T))
print('T=', isent.stat_temp_from_mach(M, Tt))
print('T=', isent.stat_temp_from_vel(v, Tt, 'True'))
print('T=', isent.stat_temp_from_vel(v, Tt, 'False'))
print('p=', isent.stat_pressure_from_mach(M, pt, Tt))

print('--------------------------------------------------------')

from pyturb.gas_models import isa


perf_air = PerfectIdealGas('Air')
perf_air_isentflow = IsentropicFlow(perf_air)

h0 = 5000 #m
T0 = isa.temperature_isa(h0)
p0 = isa.pressure_isa(h0)

rho0 = p0/perf_air_isentflow.fluid.Rg/T0

print("At {0:4.0f}m above SL:  T0={1:6.2f} K;  p0={2:6.0f} Pa;  rho0={3:6.4f} kg/m^3".format(h0, T0, p0, rho0))

mach = 0.2

vel = perf_air_isentflow.vel_from_mach(mach, T0)
T0t = perf_air_isentflow.stag_temp_from_mach(mach, T0)

qi_p = perf_air_isentflow.dynamic_pressure_from_statpress(mach, p0)
qi_rho = perf_air_isentflow.dynamic_pressure(vel, rho0)

qi_dens = 0.5*rho0*vel**2
qi_pres = 0.5*perf_air.gamma()*mach*mach*p0

print('qi static pressure and Mach: qi={}Pa'.format(qi_p))
print('qi density and velocity: qi={}Pa'.format(qi_rho))
print(qi_dens)
print(qi_pres)