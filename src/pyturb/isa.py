#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
International Standard Atmosphere conditions for temperature and pressure as a
function of the geometric altitude.

Troposphere layer implemented.

Content:
--------
    + temp_ISA
    + pres_ISA
    + height_ISA

Part of pyturb. Tests pending.

April 2017

MRodriguez
"""

import air_model
import units


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
