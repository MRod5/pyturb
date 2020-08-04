"""
pyturb
air_model.py tests

M Rodriguez
"""


from pyturb.air_model import RealGas, IdealGas


ig = IdealGas()
print(ig.cp_air())
print(ig.cv_air())
print(ig.gamma_air())

rg = RealGas(cp_option='naca', gamma_option='naca')
print(rg.cp_air(273),  rg.gamma_air(273), rg.cv_air(273))

rg = RealGas(cp_option='nasa', gamma_option='naca')
print(rg.cp_air(273),  rg.gamma_air(273), rg.cv_air(273))

rg = RealGas('nasa', 'standard')
print(rg.cp_air(2000), rg.gamma_air(2000))

rg = RealGas('low-pressure')
print(rg.cp_air(273))

rg = RealGas('aaa', 'bbb')
print(rg.cp_air(273), rg.gamma_air(273))

print(rg.R_air)