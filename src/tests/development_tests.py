import numpy as np
import pyturb.utils.constants as cts
import pyturb.utils.units as units
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.thermo_properties import ThermoProperties
print('tests_1')
print('-----------------')
#print(cts.Ru, units.m_to_nm)
#print()
#tp = ThermoProperties('O2')
#print(tp)


pig = PerfectIdealGas('Air')
print(pig.Ru)
print(pig.Rg)
