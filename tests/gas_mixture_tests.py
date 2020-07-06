from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.gas_mixture import GasMixture
import numpy as np

air = PerfectIdealGas('Air')

O2 = PerfectIdealGas('O2')
N2 = PerfectIdealGas('N2')

gm = GasMixture(gas_model='Perfect')

gm.add_gas('O2', 21)
gm.add_gas('N2', 79)

print(gm.mixture_gases)
