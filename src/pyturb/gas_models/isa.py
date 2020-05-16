"""



MRodriguez. 2020
"""

import pyturb.utils.units as units
from pyturb.air_model import IdealGas
import numpy as np
import pandas as pd
import os

isa_dir = os.path.dirname(__file__) 
file_atmos1975 = os.path.join(isa_dir, r'./atmos_layers_coesa1975.dat')

coesaRu = 8.31432 # J/mol/K (universal gas constant as in COESA 1975)
Mair = 0.0289644 # kg/mol (mean molecular mass of air as in COESA 1975)
coesaRair = coesaRu/Mair


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
        layer = 'layer not implemented for h={0}m'.format(height)

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

    # Check is height is array/list or discrete value:
    if type(height) in [np.ndarray, list]:
        # Array or list of heights:
        temperature = np.zeros_like(height)

        for ii, h in enumerate(height):
            if type(h) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
                h = np.float32(h) if (type(h)==int) else h # Avoid integers to avoid truncation
                
                # Retrieve layer information:
                temp_gradient, base_temperature, _, base_height, _ = get_atmosdata(h)
                
                #Get temperature value
                T = base_temperature + temp_gradient*(h-base_height)
                temperature[ii] = T

        else:
            # If height is not aray, list, float, int...
            print('Input height ({}) is not float or np.ndarray'.format(height))
            temperature = np.nan


    elif type(height) in [float, int, np.float64, np.float32, np.int64, np.int32]:
        # Height is discrete value:
        # Retrieve layer information:
        temp_gradient, base_temperature, _, base_height, _ = get_atmosdata(height)

        # Get temperature
        temperature = base_temperature + temp_gradient*(height-base_height)

    else:
        # If height is not aray, list, float, int...
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

    # Check is height is array/list or discrete value:
    if type(height) in [np.ndarray, list]:
        # Array or list of heights:
        isa_dev = np.zeros_like(height) if isa_dev==0 else isa_dev

        pressure = np.zeros_like(height)

        for ii, (h, isa_dev_) in enumerate(zip(height, isa_dev)):
            if type(h) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
                # Retrieve layer information:
                temp_gradient, base_temperature, base_pressure, base_height, layer = get_atmosdata(h)
                
                # Calculate temperature at h
                temperature = temperature_isa(h, isa_dev_)

                # Check if layer is isothermal:
                if layer in ['tropopause', 'stratopause']:
                    factor = np.exp(-units.grav/coesaRair/temperature*(h-base_height))

                else:
                    factor = (temperature/base_temperature)**(units.grav/coesaRair/(-temp_gradient))

                # Get pressure value
                pressure[ii] = base_pressure * factor

            else:
                # If height is not aray, list, float, int...
                print('Input height ({}) is not float or np.ndarray'.format(height))
                pressure = np.nan


    elif type(height) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
        # Height is discrete value:
        # Retrieve layer information:
        temp_gradient, base_temperature, base_pressure, base_height, layer = get_atmosdata(height)

        # Calculate temperature at height:
        temperature = temperature_isa(height, isa_dev)

        # Check if layer is isothermal:
        if layer in ['tropopause', 'stratopause']:
            factor = np.exp(-units.grav/coesaRair/temperature*(height-base_height))
            
        else:
            factor = (temperature/base_temperature)**(units.grav/coesaRair/(-temp_gradient))

        # Get pressure value
        pressure = base_pressure * factor


    else:
        # If height is not aray, list, float, int...
        print('Input height ({}) is not float or np.ndarray'.format(height))
        pressure = np.nan


    return pressure


def height_isa(p0):

    return None


def density_isa(height):

    return None