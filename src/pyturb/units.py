"""

April 2017

MRodriguez
"""


from numpy import pi

## Constants:
grav = 9.80665      # Gravitational acceleration [m/s^2]
Ru = 8.314472       # Universal gas constant [J/mol/K]
p_ref_SL = 101325   # Reference pressure at sea level [Pa]


## Units conversion: 
#Angles:
deg2rad = pi/180    # Degrees to radians
rad2deg = 180/180   # Radians to degrees


# Distance:
ft2m = 0.3048       # feet to meters
m2ft = 1/ft2m       # meters to ft
FL2ft = 100         # Flight level to ft
in2m = 0.0254       # Inches to meters
m2in = 1/in2m       # Meters to inches
nm2m = 1852         # Nautical mile to meters
m2nm = 1/nm2m       # Meters to nautical miles


# Speed:
kts2ms = 1852/3600  # A Knot is a nautical mile per hour. Conversion to meters per second is exact this way.
ms2kts = 1/kts2ms   # Meters per second to knots
mph2ms = 1609.344/3600 # Statute  mile per hour to meters per second
ms2mph = 1/mph2ms # Meters per second to statute miles per hour


# Mass:
lb2kg = 0.45359237  # Pounds (Avoirdupois) to kilograms
kg2lb = 1/lb2kg     # Kilograms to Pounds (Avoirdupois)
slug2kg = grav*m2ft*lb2kg # Slugs to kilograms
kg2slug = 1/slug2kg # Kilograms to slugs


# Force:
kp2N = grav         # Kilopond to Newton
N2kp = 1/grav       # Newtons to Kiloponds
lbf2kp = lb2kg      # Pound-force to kiloponds
kp2lb = 1/lb2kg     # Kiloponds to pound-force
lbf2N = lbf2kp*kp2N # Pound-force to Newtons
N2lbf = 1/lbf2N     # Pound-force to Newtons


# Pressure:
bar2Pa = 1e5        # Bar to Pascals
Pa2bar = 1e-5       # Pascals to bar
mbar2Pa = 1e-3*bar2Pa # mbar (mili bar) to Pascal
Pa2mbar = 1e+3*Pa2bar # Pascal to mbar (mili bar)
atm2Pa = p_ref_SL   # Atmospheres to Pascals
Pa2atm = 1/atm2Pa   # Pascals to Atmospheres
Torr2Pa = p_ref_SL/760 # Torr (mmHg) to Pascals
Pa2Torr = 1/Torr2Pa # Torr (mmHg) to Pascals
mmHg2Pa = Torr2Pa   # Torr = mmHg
Pa2mmHg = Pa2Torr   # Torr = mmHg


# Energy:
kWh2J = 3.6e6       # kilowatt-hour to Joules
J2kWh = 1/kWh2J     # Joules to kilowatt-hour
cal2J = kWh2J/860*1e-3 # Calories (international steam table) to joules
J2cal = 1/cal2J     # Joules to calories (international steam table)
Btu2J =  1055.06    # BTU (British thermal unit - ISO) to Joules
J2Btu = 1/Btu2J     # Joules to BTU (ISO)


# Power:
hp2W = 550 * lbf2N * ft2m # Horsepower (mechanical) to Watts
W2hp = 1/hp2W       # Watts to Horsepower (mechanical)
metric_hp2W = 75 * kp2N  # Horsepower (metric) to Watts
W2hp_metric = 1/metric_hp2W       # Watts to Horsepower (metric)


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
