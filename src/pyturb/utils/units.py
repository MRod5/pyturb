"""

Universal constants and Conversion of units to and from International System of units.

Content:
--------
    Constants:
    ----------
        grav : Gravitational acceleration
        Ru : Universal gas constant
    
    Conversion of units:
    --------------------
        Angles:
            Radians and degrees

        Distance: Between meters [m] and:
            feet (ft), inches (in), Flight Level (FL), Nautical Mile (nm)
        
        Speed: Between meters per secons (ms) [m/s] and:
            Knots [kts], mile sper hour [MPH]
        
        Mass: Between kilograms [kg] and:
            Pounds [lb], slugs [slug]
        
        Force: Between Newtons [N] and:
            kiloponds [kp, kgf], pound-force [lbf]
        
        Pressure: Between Pascals [Pa] and:
            Bars [bar], milibars [mbar], Torricelli [Torr, mmHg], atmospheres [atm]
        
        Energy: Between Joules [J] and:
            kilowatt-hour [kWh], calories [cal], Briish Thermal Units [Btu]
        
        Power: Between Watts [W] and:
            metric horsepower [hp], mechanical horsepower [hp]
        
        Temperatura: Between Celsius [ºC], Kelvin [K], Rankine [R] and Fahrenheit [ºF]


MRodriguez. 2020


"""


from numpy import pi
import pyturb.utils.constants as cts


## Units conversion: 
#Angles:
deg_to_rad = pi/180     # Degrees to radians
rad_to_deg = 180/180    # Radians to degrees


# Distance:
ft_to_m = 0.3048        # feet to meters
m_to_ft = 1/ft_to_m     # meters to ft
FL_to_ft = 100          # Flight level to ft
FL_to_m = FL_to_ft*ft_to_m # Flight level to ft
in_to_m = 0.0254        # Inches to meters
m_to_in = 1/in_to_m     # Meters to inches
nm_to_m = 1852          # Nautical mile to meters
m_to_nm = 1/nm_to_m     # Meters to nautical miles


# Speed:
kts_to_ms = 1852/3600   # A Knot is a nautical mile per hour. Conversion to meters per second is exact this way.
ms_to_kts = 1/kts_to_ms   # Meters per second to knots
mph_to_ms = 1609.344/3600 # Statute  mile per hour to meters per second
ms_to_mph = 1/mph_to_ms # Meters per second to statute miles per hour


# Mass:
lb_to_kg = 0.45359237   # Pounds (Avoirdupois) to kilograms
kg_to_lb = 1/lb_to_kg   # Kilograms to Pounds (Avoirdupois)
slug_to_kg = cts.grav*m_to_ft*lb_to_kg # Slugs to kilograms
kg_to_slug = 1/slug_to_kg # Kilograms to slugs


# Force:
kp_to_N = cts.grav          # Kilopond to Newton
N_to_kp = 1/cts.grav        # Newtons to Kiloponds
lbf_to_kp = lb_to_kg    # Pound-force to kiloponds
kp_to_lb = 1/lb_to_kg   # Kiloponds to pound-force
lbf_to_N = lbf_to_kp*kp_to_N # Pound-force to Newtons
N_to_lbf = 1/lbf_to_N   # Pound-force to Newtons


# Pressure:
bar_to_Pa = 1e5         # Bar to Pascals
Pa_to_bar = 1e-5        # Pascals to bar
mbar_to_Pa = 1e-3*bar_to_Pa # mbar (mili bar) to Pascal
Pa_to_mbar = 1e+3*Pa_to_bar # Pascal to mbar (mili bar)
atm_to_Pa = cts.p_ref_SL    # Atmospheres to Pascals
Pa_to_atm = 1/atm_to_Pa # Pascals to Atmospheres
Torr_to_Pa = cts.p_ref_SL/760 # Torr (mmHg) to Pascals
Pa_to_Torr = 1/Torr_to_Pa # Torr (mmHg) to Pascals
mmHg_to_Pa = Torr_to_Pa # Torr = mmHg
Pa_to_mmHg = Pa_to_Torr # Torr = mmHg


# Energy:
kWh_to_J = 3.6e6        # kilowatt-hour to Joules
J_to_kWh = 1/kWh_to_J   # Joules to kilowatt-hour
cal_to_J = kWh_to_J/860*1e-3 # Calories (international steam table) to joules
J_to_cal = 1/cal_to_J   # Joules to calories (international steam table)
Btu_to_J =  1055.06     # BTU (British thermal unit - ISO) to Joules
J_to_Btu = 1/Btu_to_J   # Joules to BTU (ISO)


# Power:
hp_to_W = 550 * lbf_to_N * ft_to_m  # Horsepower (mechanical) to Watts
W_to_hp = 1/hp_to_W                 # Watts to Horsepower (mechanical)
metric_hp_to_W = 75 * kp_to_N       # Horsepower (metric) to Watts
W_to_hp_metric = 1/metric_hp_to_W   # Watts to Horsepower (metric)


# Temperature:
def celsius_to_kelvin(temp):
    """
    From Celsius (ºC) to Kelvin (K)
    """
    return temp + 273.15


def kelvin_to_celsius(temp):
    """
    From Kelvin (K) to Celsius (ºC)
    """
    return temp - 273.15


def celsius_to_fahrenheit(temp):
    """
    From Celsius (ºC) to Fahrenheit (ºF)
    """
    return temp*9/5 + 32


def fahrenheit_to_celsius(temp):
    """
    From Fahrenheit (ºF) to Celsius (ºC)
    """
    return (temp - 32)*5/9


def fahrenheit_to_rankine(temp):
    """
    From Fahrenheit (ºF) to Rankine (R)
    """
    return temp + 459.67


def rankine_to_fahrenheit(temp):
    """
    From Rankine (R) to Fahrenheit (ºF)
    """
    return temp - 459.67


def fahrenheit_to_kelvin(temp):
    """
    From Fahrenheit (ºF) to Kelvin (K)
    """
    return celsius_to_kelvin(fahrenheit_to_celsius(temp))


def kelvin_to_fahrenheit(temp):
    """
    From Kelvin (K) to Fahrenheit (ºF)
    """
    return celsius_to_fahrenheit(kelvin_to_celsius(temp))


def celsius_to_rankine(temp):
    """
    From Celsius (ºC) to Rankine (R)
    """
    return fahrenheit_to_rankine(celsius_to_fahrenheit(temp))


def rankine_to_celsius(temp):
    """
    From Rankine (R) to Celsius (ºC)
    """
    return fahrenheit_to_celsius(rankine_to_fahrenheit(temp))


def kelvin_to_rankine(temp):
    """
    From Kelvin (K) to Rankine (R)
    """
    return temp*9/5


def rankine_to_kelvin(temp):
    """
    From Rankine (R) to Kelvin (K)
    """
    return temp*5/9

