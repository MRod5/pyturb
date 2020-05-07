#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities and auxiliary functions

Content:
--------
    
Dependencies:
-------------


Part of pyturb. Tests pending.

April 2017

MRodriguez
"""

# *** Constants:
grav = 9.80665  # Acceleration of gravity [m/s^2]

# *** Units conversion:

# Distance:
ft2m = 0.3048  # foot to meter
m2ft = 1/ft2m  # meters to ft
kts2ms = 1852/3600  # A Knot is a mile per hour. Conversion to meters per second is exact this way.
ms2kts = 1/kts2ms  # meters per second to knots


# Temperature:
def celsius2kelvin(temp):
    return temp + 273.15


def kelvin2celsius(temp):
    return temp - 273.15


def celsius2fahrenheit(temp):
    return temp*9/5 + 32


def fahrenheit2celsius(temp):
    return (temp - 32)*5/9


def fahrenheit2rankine(temp):
    return temp + 459.67


def rankine2fahrenheit(temp):
    return temp - 459.67


def fahrenheit2kelvin(temp):
    return celsius2kelvin(fahrenheit2celsius(temp))


def kelvin2fahrenheit(temp):
    return celsius2fahrenheit(kelvin2celsius(temp))


def celsius2rankine(temp):
    return fahrenheit2rankine(celsius2fahrenheit(temp))


def rankine2celsius(temp):
    return fahrenheit2celsius(rankine2fahrenheit(temp))


def kelvin2rankine(temp):
    return temp*9/5


def rankine2kelvin(temp):
    return temp*5/9
