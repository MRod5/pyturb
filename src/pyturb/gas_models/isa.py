"""
ISA: International Standard Atmosphere:
---------------------------------------

The International Standard Atmosphere stablishes a model for the temperature, 
pressure, density and Mach number as a function of the geopotential heights.

The ISA model refers to the US Atmosphere model of 1975 (COESA 1975), although 
from 50 km on, ISA model stablishes the 1975 US Atmosphere Model as interim.

Content:
--------
    + get_atmosdata: provides atmosphere layer information given a geopotential heights
    + temperature_isa: ISA temperature as a function of the geopotential height
    + pressure_isa: ISA pressure as a function of the geopotential height
    + density_isa: ISA density as a function of the geopotential height
    + density_state_eq: Density value from the ideal gas law. p and T are calculated 
                        with the ISA coefficients gor a given geopotential height
    + height_from_temperature_isa: geopotential height from the temperature assuming
                                   ISA coefficients.
    + height_from_pressure_isa: geopotential height from the pressure assuming ISA coefficients

References:
-----------
- "Defining constants, equations and abbreviated tables of the 1975 U.S. Standard Atmosphere" NASA TR R-459
- "Computational modelling of aircraft and the envirorenment" Dominic J. Diston


MRodriguez. 2020

"""

import pyturb.utils.units as units
import pyturb.utils.constants as cts
import numpy as np
import pandas as pd
import os


# Directory to the 1975 COESA data:
isa_dir = os.path.dirname(__file__) 
file_atmos1975 = os.path.join(isa_dir, r'./atmos_layers_coesa1975.dat')

# Ideal gas constant and molecular mass of air as in 1975 COESA
coesaRair = cts.coesaRu/cts.coesaMair


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
        base_density = atmos_data.iloc[layer_mask[0][-1]]['base_density']
        base_height = atmos_data.iloc[layer_mask[0][-1]]['geopotential_height']
        layer = atmos_data.iloc[layer_mask[0][-1]]['atmos_layer'].strip()
    else:
        # Layer not implemented:
        temp_gradient = np.nan
        base_temperature = np.nan
        base_pressure = np.nan
        base_density = np.nan
        base_height = np.nan
        layer = 'layer not implemented for h={0}m'.format(height)

    return temp_gradient, base_temperature, base_pressure, base_density, base_height, layer


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
                 isa_dev can be:
                 + A discrete value: the same isa_dev is applied all temperatures
                 + An array: every temperature has its own deviation
                 By default isa_dev=0
        
    + Outputs:
    ----------
        temperature: ndarray or float. Static temperature at input height [K]
        
    """

    # Check if height is array/list or discrete value:
    if type(height) in [np.ndarray, list]:
        # Array or list of heights:
        temperature = np.zeros_like(height, dtype=np.float64)

        isa_dev = isa_dev*np.ones_like(height, dtype=np.float64) if np.size(isa_dev)==1 else isa_dev
        for ii, (h, isa_dev_) in enumerate(zip(height, isa_dev)):
            if type(h) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
                # Retrieve layer information:
                temp_gradient, base_temperature, _, _, base_height, _ = get_atmosdata(h)
                
                #Get temperature value
                T = base_temperature + temp_gradient*(h-base_height)
                temperature[ii] = T + isa_dev_

            else:
                # If height is not aray, list, float, int...
                print('Input height ({}) is not numeric or np.ndarray'.format(height))
                temperature = np.nan


    elif type(height) in [float, int, np.float64, np.float32, np.int64, np.int32]:
        # Height is discrete value:
        # Retrieve layer information:
        temp_gradient, base_temperature, _, _, base_height, _ = get_atmosdata(height)

        # Get temperature
        temperature = base_temperature + temp_gradient*(height-base_height)
        temperature += isa_dev

    else:
        # If height is not aray, list, float, int...
        print('Input height ({}) is not numeric or np.ndarray'.format(height))
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
        isa_dev: float. Should be zero for pressure-altitude with ISA coefficients
                 Standard day base temperature deviation [K]
                 isa_dev can be:
                 + A discrete value: the same isa_dev is applied all temperatures
                 + An array: every temperature has its own deviation
                 By default isa_dev=0
        
    + Outputs:
    ----------
        pressure: ndarray or float. Static pressure at input height [Pa]

    """

    # Check if height is array/list or discrete value:
    if type(height) in [np.ndarray, list]:
        # Array or list of heights:
        isa_dev = isa_dev*np.ones_like(height, dtype=np.float64) if np.size(isa_dev)==1 else isa_dev

        pressure = np.zeros_like(height, dtype=np.float64)

        for ii, (h, isa_dev_) in enumerate(zip(height, isa_dev)):
            if type(h) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
                # Retrieve layer information:
                temp_gradient, base_temperature, base_pressure, _, base_height, layer = get_atmosdata(h)
                
                # Calculate temperature at h
                temperature = temperature_isa(h, isa_dev_)

                # Check if layer is isothermal:
                if layer in ['tropopause', 'stratopause', 'mesosphere3']:
                    # XXX Posible errata en la temperatura, debe ser temp_base
                    factor = np.exp(-cts.grav/coesaRair/temperature*(h-base_height))

                else:
                    factor = (temperature/base_temperature)**(cts.grav/coesaRair/(-temp_gradient))

                # Get pressure value
                pressure[ii] = base_pressure * factor

            else:
                # If height is not aray, list, float, int...
                print('Input height ({}) is not numeric or np.ndarray'.format(height))
                pressure = np.nan


    elif type(height) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
        # Height is discrete value:
        # Retrieve layer information:
        temp_gradient, base_temperature, base_pressure, _, base_height, layer = get_atmosdata(height)

        # Calculate temperature at height:
        temperature = temperature_isa(height, isa_dev)

        # Check if layer is isothermal:
        if layer in ['tropopause', 'stratopause', 'mesosphere3']:
            factor = np.exp(-cts.grav/coesaRair/temperature*(height-base_height))
            
        else:
            factor = (temperature/base_temperature)**(cts.grav/coesaRair/(-temp_gradient))

        # Get pressure value
        pressure = base_pressure * factor


    else:
        # If height is not aray, list, float, int...
        print('Input height ({}) is not numeric or np.ndarray'.format(height))
        pressure = np.nan


    return pressure


def density_isa(height, isa_dev=0):
    """
    ISA density:
    -------------
    
    International Standard Atmosphere density-altitude as a function of
    the geopotential height.

    + Inputs:
    ---------
        h: ndarray or list or float or int. Geometric altitude [m]
        isa_dev: float. Standard day base temperature deviation [K]

        
    + Outputs:
    ----------
        density: ndarray or float. Static density at input height [Pa]

    """

    # Check if height is array/list or discrete value:
    if type(height) in [np.ndarray, list]:
        # Array or list of heights:
        isa_dev = isa_dev*np.ones_like(height, dtype=np.float64) if np.size(isa_dev)==1 else isa_dev

        density = np.zeros_like(height, dtype=np.float64)

        for ii, (h, isa_dev_) in enumerate(zip(height, isa_dev)):
            if type(h) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
                # Retrieve layer information:
                temp_gradient, base_temperature, _, base_density, base_height, layer = get_atmosdata(h)
                
                # Calculate temperature at h
                temperature = temperature_isa(h, isa_dev_)

                # Check if layer is isothermal:
                if layer in ['tropopause', 'stratopause', 'mesosphere3']:
                    factor = np.exp(-cts.grav/coesaRair/temperature*(h-base_height))

                else:
                    factor = (temperature/base_temperature)**(cts.grav/coesaRair/(-temp_gradient))

                # Get density value
                density[ii] = base_density * factor

            else:
                # If height is not aray, list, float, int...
                print('Input height ({}) is not numeric or np.ndarray'.format(height))
                density = np.nan


    elif type(height) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
        # Height is discrete value:
        # Retrieve layer information:
        temp_gradient, base_temperature, _, base_density, base_height, layer = get_atmosdata(height)

        # Calculate temperature at height:
        temperature = temperature_isa(height, isa_dev)

        # Check if layer is isothermal:
        if layer in ['tropopause', 'stratopause', 'mesosphere3']:
            factor = np.exp(-cts.grav/coesaRair/temperature*(height-base_height))
            
        else:
            factor = (temperature/base_temperature)**(cts.grav/coesaRair/(-temp_gradient))

        # Get density value
        density = base_density * factor


    else:
        # If height is not aray, list, float, int...
        print('Input height ({}) is not numeric or np.ndarray'.format(height))
        density = np.nan


    return density


def density_state_eq(height, isa_dev=0):
    """
    ISA density:
    -------------
    
    International Standard Atmosphere density-altitude from the ideal gas law equation:
        p = rho * Rair * T

    + Inputs:
    ---------
        h: ndarray or list or float or int. Geometric altitude [m]
        isa_dev: float. Standard day base temperature deviation [K]
        
    + Outputs:
    ----------
        density: ndarray or float. Static density at input height [Pa]

    """

    # Check if height is array/list or discrete value:
    if type(height) in [np.ndarray, list]:
        # Array or list of heights:
        density = np.zeros_like(height, dtype=np.float64)
        
        isa_dev = isa_dev*np.ones_like(height, dtype=np.float64) if np.size(isa_dev)==1 else isa_dev

        for ii, (h, isa_dev_) in enumerate(zip(height, isa_dev)):
            if type(h) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
                p = pressure_isa(h, isa_dev_)
                T = temperature_isa(h, isa_dev_)
                rho = p/coesaRair/T

                density[ii] = rho

    elif type(height) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
        # Height is discrete value:
                p = pressure_isa(height, isa_dev)
                T = temperature_isa(height, isa_dev)
                rho = p/coesaRair/T

                density = rho

    else:
        # If height is not aray, list, float, int...
        print('Input height ({}) is not float or np.ndarray'.format(height))
        density = np.nan


    return density


def height_from_temperature_isa(temperature, isa_dev=0, layer = None):
    """
    height_from_temperature_isa:
    ----------------------------

    Geopotential height that corresponds the the input temperature value, regarding the International
    Standard Atmosphere.

    Note that for a given temperature, depending on the layer, more than one geopotential height can be
    obtained. The output array, therefore, holds all the heights that can be obtained for a given temperature.
    Thus, note that if the temperature value produces more than one geopotential height, an array of arrays 
    with be returned, with as many rows as temperatures provided and as many columns as heights compatible 
    with the temperature value.

    + Inputs:
    ---------
        temperature: ndarray or list or float or int. Temperature in [k]
        isa_dev: float. Standard day base temperature deviation [K]
        
    + Outputs:
    ----------
        height: ndarray or float. Geopotential heights corresponding to the input temperature [m]

    """

    # Check if temperature is array/list or discrete value:
    if type(temperature) in [np.ndarray, list]:
        # Array or list of temperatures:
        isa_dev = isa_dev*np.ones_like(height, dtype=np.float64) if np.size(isa_dev)==1 else isa_dev
        
        height = []

        for T, isa_dev_ in zip(temperature, isa_dev):
            if type(T) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
                atmos_data = pd.read_csv(file_atmos1975, sep='|')

                # Temperature without isa deviation
                T_isa = T - isa_dev_

                temp_gradient = atmos_data['temperature_gradient']
                base_temperature = atmos_data['base_temperature']
                base_height = atmos_data['geopotential_height']
                height_limit = atmos_data['geopotential_height_limit']
                layer = atmos_data['atmos_layer']
        
                # For each layer the geopotential height is calculated
                height_ = []
                for temp_grad, base_temp, base_h, layer_, hlim in zip(temp_gradient, base_temperature, base_height, layer, height_limit):
                    if layer_.strip() not in ['tropopause', 'stratopause', 'mesosphere3']:
                        # Height value:
                        h = (T_isa - base_temp)/temp_grad + base_h

                        # If height value is compatible with the layer store it
                        if base_h<=h<=hlim:
                            height_.append(h)
                
                # Output heights for current temperature:
                height_ = np.asarray(height_[:])
                height.append(height_)
        
        # Output height array:
        height = np.asarray(height)
        
    elif type(temperature) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
         # Temperature is a discrete value:

        if layer is None:
            atmos_data = pd.read_csv(file_atmos1975, sep='|')

            temperature_isa = temperature - isa_dev

            temp_gradient = atmos_data['temperature_gradient']
            base_temperature = atmos_data['base_temperature']
            base_height = atmos_data['geopotential_height']
            height_limit = atmos_data['geopotential_height_limit']
            layer = atmos_data['atmos_layer']

            height = []
            for temp_grad, base_temp, base_h, layer_, hlim in zip(temp_gradient, base_temperature, base_height, layer, height_limit):
                if layer_.strip() not in ['tropopause', 'stratopause', 'mesosphere3']:
                    # Height value:
                    h = (temperature_isa - base_temp)/temp_grad + base_h
                    
                    # If height value is compatible with the layer store it
                    if base_h<=h<=hlim:
                        height.append(h)

            # Output heights for current temperature:
            height = np.asarray(height)

        else:
            # If layer is known (obtained from pressure), only one height value can be computed
            atmos_data = pd.read_csv(file_atmos1975, sep='|')

            # Filtered layer data:
            base_temp = atmos_data.loc[atmos_data['atmos_layer']==layer]['base_temperature'].values[0]
            temp_grad = atmos_data.loc[atmos_data['atmos_layer']==layer]['temperature_gradient'].values[0]
            base_h = atmos_data.loc[atmos_data['atmos_layer']==layer]['geopotential_height'].values[0]
            
            # Height value
            temperature_isa = temperature - isa_dev
            height = (temperature_isa - base_temp)/temp_grad + base_h
        

    else:
        # If pressure is not aray, list, float, int...
        print('Input pressure ({}) is not numeric or np.ndarray'.format(height))
        height = np.nan
    
    return height


def height_from_pressure_isa(pressure, isa_dev=0):
    """
    height_from_pressue_isa:
    ------------------------

    Geopotential height that corresponds the the input pressure value, regarding the International
    Standard Atmosphere.

    + Inputs:
    ---------
        pressure: ndarray or list or float or int. Altitude-pressure [Pa]
        isa_dev: float. Standard day base temperature deviation [K]
        
    + Outputs:
    ----------
        height: ndarray or float. Geopotential height corresponding to the input pressure [m]

    """

    # Check if pressure is array/list or discrete value:
    if type(pressure) in [np.ndarray, list]:
        # Array or list of pressure:
        isa_dev = np.zeros_like(pressure, dtype=np.float64) if isa_dev==0 else isa_dev
        
        height = np.zeros_like(pressure, dtype=np.float64)

        for ii, (p, isa_dev_) in enumerate(zip(pressure, isa_dev)):
            if type(p) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
                # Retrieve layer information:
                atmos_data = pd.read_csv(file_atmos1975, sep='|')
                layer_mask = np.where(atmos_data['base_pressure'].values>=p)
                
                # Layer mask where pressure is located:
                temp_gradient = atmos_data.iloc[layer_mask[0][-1]]['temperature_gradient']
                base_temperature = atmos_data.iloc[layer_mask[0][-1]]['base_temperature']
                base_pressure = atmos_data.iloc[layer_mask[0][-1]]['base_pressure']
                base_height = atmos_data.iloc[layer_mask[0][-1]]['geopotential_height']
                layer = atmos_data.iloc[layer_mask[0][-1]]['atmos_layer']
        
            if layer.strip() in ['tropopause', 'stratopause', 'mesosphere3']:
                h = base_height - np.log(p/base_pressure)*coesaRair/cts.grav*base_temperature

            else:
                Tratio = (p/base_pressure)**(-temp_gradient*coesaRair/cts.grav)
                T = Tratio * base_temperature
                h = height_from_temperature_isa(T, isa_dev_, layer)
        
            height[ii] = h

    elif type(pressure) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
         # Pressure is a discrete value:

        if type(pressure) in  [float, int, np.float64, np.float32, np.int64, np.int32]:
            # Retrieve layer information:
            atmos_data = pd.read_csv(file_atmos1975, sep='|')
            layer_mask = np.where(atmos_data['base_pressure'].values>=pressure)
                
            # Layer mask where pressure is located:
            temp_gradient = atmos_data.iloc[layer_mask[0][-1]]['temperature_gradient']
            base_temperature = atmos_data.iloc[layer_mask[0][-1]]['base_temperature']
            base_pressure = atmos_data.iloc[layer_mask[0][-1]]['base_pressure']
            base_height = atmos_data.iloc[layer_mask[0][-1]]['geopotential_height']
            layer = atmos_data.iloc[layer_mask[0][-1]]['atmos_layer']
        
        if layer.strip() in ['tropopause', 'stratopause', 'mesosphere3']:
            # If layer is isothermal
            h = base_height - np.log(pressure/base_pressure)*coesaRair/cts.grav*base_temperature

        else:
            # If temperature gradient is not zero:
            Tratio = (pressure/base_pressure)**(-temp_gradient*coesaRair/cts.grav)
            T = Tratio * base_temperature
            h = height_from_temperature_isa(T, isa_dev, layer)
        
        height = h

    else:
        # If pressure is not aray, list, float, int...
        print('Input pressure ({}) is not numeric or np.ndarray'.format(height))
        height = np.nan

    return height