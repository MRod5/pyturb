from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.gas_mixture import GasMixture
import numpy as np

air = PerfectIdealGas('Air')

O2 = PerfectIdealGas('O2')
N2 = PerfectIdealGas('N2')

gm = GasMixture(gas_model='Perfect')

gm.add_gas('O2', 0.209476)
gm.add_gas('N2', 0.78084)
gm.add_gas('Ar', 0.009365)
gm.add_gas('CO2', 0.000319)

print(gm.mixture_gases)
print('-------------------------------------------------------')
print('Ng',gm.mixture_gases['Ng'], np.sum(gm.mixture_gases['Ng']))
print('Mg',gm.Mg, 'Rg', gm.Rg)

print(gm.mixture_gases.loc[0]['gas_species'], gm.mixture_gases.loc[0]['gas_properties'].cp(288))
print(gm.mixture_gases.loc[1]['gas_species'], gm.mixture_gases.loc[1]['gas_properties'].cp(288))
print('-------------------------------------------------------')

np.testing.assert_almost_equal(gm.Rg, air.Rg)
np.testing.assert_almost_equal(gm.Mg, air.thermo_prop.Mg)