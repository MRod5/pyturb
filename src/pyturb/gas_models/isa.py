"""



MRodriguez. 2020
"""

import pyturb.air_model as air_model
import pyturb.utils.units as units
import numpy as np
import pandas as pd

file_atmos1975 = './atmos_layers_coesa1975.dat'


def get_atmosdata(h):
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
    atmos_data = pd.read_csv(file_atmos1975, sep='\t')
    layer_mask = np.where(atmos_data['geopotential_height'].values>h)

    if np.size(layer_mask)>0:
        # Extract infomation of the layer corresponding to h:
        temp_gradient = atmos_data.iloc[layer_mask[0][0]]['temperature_gradient']
        base_temperature = atmos_data.iloc[layer_mask[0][0]]['base_temperature']
        base_pressure = atmos_data.iloc[layer_mask[0][0]]['base_pressure']
        base_height = atmos_data.iloc[layer_mask[0][0]]['geopotential_height']
        layer = atmos_data.iloc[layer_mask[0][0]]['atmos_layer']
    else:
        # Layer not implemented:
        temp_gradient = np.nan
        base_temperature = np.nan
        base_pressure = np.nan
        base_height = np.nan
        layer = 'layer not implemented for h={0}m'.format(h)

    return temp_gradient, base_temperature, base_pressure, base_height, layer


def temp_isa(h, isa_dev=0):
    """
    ISA temperature:
    ----------------
    
    International Standard Atmosphere temperature as a function of the
    altitude above sea level.
    
    + Inputs:
    ---------
        h: float. Geometric altitude [m]
        isa_dev: float. Standard day base temperature deviation [K]
        
    + Outputs:
    ----------
        temp_isa: float. Static temperature at altitude "h" [K]
        
    """
    if h <= 11000:
        
        temp_base = 15 + 273.15 + isa_dev  # Base temperature at sea level [K]
        temp_rate = 6.5e-3  # Layer rate of temperature [K/m]

        return temp_base - temp_rate*h
    else:
        print("Layer not implemented yet for altitude {}".format(h))
        return
    

def pres_isa(h):
    """
    ISA pressure:
    -------------
    
    International Standard Atmosphere pressure as a function of the altitude
    above sea level.

    + Inputs:
    ---------
        h: float. Geometric altitude [m]
        
    + Outputs:
    ----------
        press_isa: float. Static pressure at altitude "h" [Pa]

    """
    if h <= 11000:
        press_base = 101325  # SL base pressure [Pa], standard day
        temp_rate = 6.5e-3  # Layer rate of temperature [K/m]
        temp_base = 15 + 273.15  # Base temperature at sea level [K], standard day

        # Get constants:
        g = units.grav
        am = air_model.Air()
        air_constant = am.R_air

        # Calculate pressure:
        press_isa = press_base*((1-temp_rate/temp_base*h)**(g/air_constant/temp_rate))
        
        return press_isa
    else:
        print("Layer not implemented yet for altitude {}".format(h))
        return


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
