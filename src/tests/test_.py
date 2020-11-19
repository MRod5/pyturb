from pyturb.gas_models.isa import temperature_isa
import pyturb.gas_models.isa as isa

T = isa.temperature_isa(1000)
T2 = temperature_isa(1000)
print(T, T2)