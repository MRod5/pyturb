from pyturb.gas_models.gas_mixture import GasMixture
from pyturb.combustion.combustion_thermodynamics import Combustion
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas

air = PerfectIdealGas('Air')
fuel = PerfectIdealGas('C8H18,isooctane')

print(fuel.thermo_prop.chemical_formula)
print(air.thermo_prop.chemical_formula)

comb = Combustion(fuel, air)

comb.combustion_stoichiometry()
print(comb.stoichiometric_reaction)

print(comb.stoich_far)
print(comb.alpha)
print(comb.beta)
print(comb.gamma)

comb.heat_of_combustion()
print(comb.hcomb_g, comb.hcomb_l)
print(comb.LHV, comb.HHV)
print(comb.products_dictionary)


oxidmix = GasMixture(gas_model="perfect")
oxidmix.add_gas('O2', 0.5)
oxidmix.add_gas('O3', 0.33333)
oxidmix.add_gas('N2', 1)

comb = Combustion(fuel, oxidmix)

comb.combustion_stoichiometry()
print(comb.stoichiometric_reaction)

comb.heat_of_combustion()
print(comb.hcomb_g, comb.hcomb_l)
print(comb.LHV, comb.HHV)
print(comb.products_dictionary)