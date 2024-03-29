"""
pyturb
Combustion thermodynamics tests

M Rodriguez. 2020
"""
import sys
from sys import path
from os.path import dirname as dir
sys.path.append(dir(sys.path[0]))

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


C2H6O = PerfectIdealGas('CH3OCH3')
O2 = PerfectIdealGas('O2')
C2H6O_comb = Combustion(C2H6O, O2)

C2H6O_comb.combustion_stoichiometry()
C2H6O_comb.heat_of_combustion()
print(C2H6O_comb.stoichiometric_reaction)
print(C2H6O_comb.hcomb_l)