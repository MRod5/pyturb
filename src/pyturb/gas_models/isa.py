"""



MRodriguez. 2020
"""

import pyturb.air_model as air_model
import pyturb.utils.units as units
from pyturb.air_model import IdealGas
import numpy as np
import pandas as pd

file_atmos1975 = './atmos_layers_coesa1975.dat'
air = IdealGas()


def get_atmosdata(height):
    """
    get_atmosdata:
    --------------

    Reads data in atmos_layes_coesa1975.dat and stores it into a Pandas dataframe.

    Data corresponds to the Nasa Technical Report "Defining constants, equations, and abbreviated 
    tables of the 1975 U.S. Standard Atmosphere".

    The dataframe stores the temperature gradient, temperature base value and pressure base value
    for each laye from sea-level (0 m) to 85 km (mesosphere).

    + Inputs:
    ---------
        h: float. Geopotential altitude [m]

    + Outputs:
    ----------
        temp_gradient: float. Temperature gradient for a given atmosphere layer [K/m]
        base_temperature: float. Temperature at the beginning of the layer [K]
        base_pressure: float. Pressure at the beginning of the layer [Pa]
        base_height: float. Geopotential height at the beginning of the layer [m]
        layer: str. Atmosphere layer
    """

    # Read atmos data:
    atmos_data = pd.read_csv(file_atmos1975, sep='|')
    layer_mask = np.where(atmos_data['geopotential_height'].values<=height)

    if np.size(layer_mask)>0:
        # Extract infomation of the layer corresponding to h:
        temp_gradient = atmos_data.iloc[layer_mask[0][-1]]['temperature_gradient']
        base_temperature = atmos_data.iloc[layer_mask[0][-1]]['base_temperature']
        base_pressure = atmos_data.iloc[layer_mask[0][-1]]['base_pressure']
        base_height = atmos_data.iloc[layer_mask[0][-1]]['geopotential_height']
        layer = atmos_data.iloc[layer_mask[0][-1]]['atmos_layer'].strip()
    else:
        # Layer not implemented:
        temp_gradient = np.nan
        base_temperature = np.nan
        base_pressure = np.nan
        base_height = np.nan
        layer = 'layer not implemented for h={0}m'.format(h)

    return temp_gradient, base_temperature, base_pressure, base_height, layer


def temperature_isa(height, isa_dev=0):
    """
    ISA temperature:
    ----------------
    
    International Standard Atmosphere temperature as a function of the
    geopotential height [m].
    
    + Inputs:
    ---------
        h: ndarray or list or float or int. Geometric altitude [m]
        isa_dev: float. Standard day base temperature deviation [K]
        
    + Outputs:
    ----------
        temperature: ndarray or float. Static temperature at input height [K]
        
    """

    if type(height) in [np.ndarray, list]:
        temperature = np.zeros_like(height)
        for ii, h in enumerate(height):
            temp_gradient, base_temperature, base_pressure, base_height, _ = get_atmosdata(h)
            T = base_temperature + temp_gradient*(h-base_height)
            temperature[ii] = T

    elif type(height) in [float, int, np.float64, np.float32, np.int64, np.int32]:
        temp_gradient, base_temperature, base_pressure, base_height, _ = get_atmosdata(height)
        temperature = base_temperature + temp_gradient*(height-base_height)

    else:
        print('Input height ({}) is not float or np.ndarray'.format(height))
        temperature = np.nan

    
    return temperature
    

def pressure_isa(height, isa_dev=0):
    """
    ISA pressure:
    -------------
    
    International Standard Atmosphere pressure-altitude as a function of
    the geopotential height.

    + Inputs:
    ---------
        h: ndarray or list or float or int. Geometric altitude [m]
        isa_dev: float. Standard day base temperature deviation [K]
        
    + Outputs:
    ----------
        pressure: ndarray or float. Static pressure at input height [Pa]

    """

    if type(height) in [np.ndarray, list]:
        pressure = np.zeros_like(height)
        if isa_dev==0:
            isa_dev = np.zeros_like(height)

        
        for ii, (h, isa_dev_) in enumerate(zip(height, isa_dev)):
            temp_gradient, base_temperature, base_pressure, base_height, layer = get_atmosdata(h)
            
            temperature = temperature_isa(h, isa_dev_)

            if layer in ['tropopause', 'stratopause']:
                factor = np.exp(-units.grav/air.R_air/temperature*(h-base_height))
                pressure[ii] = base_pressure * factor
            else:
                factor = (temperature/base_temperature)**(units.grav/air.R_air/(-temp_gradient))
                pressure[ii] = base_pressure * factor


    elif type(height) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
        temp_gradient, base_temperature, base_pressure, base_height, layer = get_atmosdata(height)

        temperature = temperature_isa(height, isa_dev)

        if layer in ['tropopause', 'stratopause']:
            factor = np.exp(-units.grav/air.R_air/temperature*(height-base_height))
            pressure = base_pressure * factor
        else:
            
            factor = (temperature/base_temperature)**(units.grav/air.R_air/(-temp_gradient))
            pressure = base_pressure * factor


    else:
        print('Input height ({}) is not float or np.ndarray'.format(height))
        pressure = np.nan

    return pressure


def height_isa(p0):
    """
    ISA Height:
    -----------
    
    Height above sea level, assuming the International Standard Atmosphere pressure model.
    
    + Inputs:
    ---------
        P0: float. Static pressure at a given altitude [Pa]
        
    + Outputs:
    ----------
        h0: float. Height over sea level [m]
    """

    press_base = 101325  # SL base pressure [Pa], standard day
    temp_rate = 6.5e-3  # Layer rate of temperature [K/m]
    temp_base = 15 + 273.15  # Base temperature at sea level [K], standard day

    # Get constants:
    g = units.grav
    am = air_model.Air()
    air_constant = am.R_air

    h0 = temp_base/temp_rate*(1 - ((p0/press_base)**(air_constant*temp_rate/g)))

    return h0


def density_isa(height):

    return None