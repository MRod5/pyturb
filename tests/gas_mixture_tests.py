from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.gas_mixture import GasMixture

O2 = PerfectIdealGas('O2')
print(O2.thermo_prop.chemical_formula)
N2 = PerfectIdealGas('N2')
print(N2.thermo_prop.chemical_formula)

gm = GasMixture(gas_model='Perfect')

gm.add_gas('O2', 1)

print(gm.mixture_gases.iloc[0])