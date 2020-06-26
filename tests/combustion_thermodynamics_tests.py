from pyturb.combustion.combustion_thermodynamics import ThermoCombustion
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas

air = PerfectIdealGas('Air')
fuel = PerfectIdealGas('C8H18,isooctane')

print(fuel.thermo_prop)

comb = ThermoCombustion(fuel, air)

comb.stoichiometry()

print(comb.stoichiometric_reaction)

