"""
Constants:
----------

Universal constants and reference values:
    + Gravity
    + Pressure at sea level
    + Temperature at sea level
    + Reference temperature for chemical reactions
    + Avogadro constant
    + Boltzmann constant
    + Universal gas law constant


MRodriguez. 2020


"""

# Gravitational acceleration
grav = 9.80665          #  [m/s^2]

# Reference pressure at sea level
p_ref_SL = 101325       # [Pa]

# Reference temperature at sea lebel
T_ref_SL = 288.15       # [K]

# Reference temperature for chemical reactions
T_ref    = 298.15       # [K]

# Avogadro constant
Na = 6.02214076e23      # [mol-1]

# Boltzmann constant
kb = 1.380649e-23       # [J/K]

# Universal gas constant
Ru = Na*kb              # 8.314462 [J/mol/K]

# Universal gas constants as in COESA
# (US Atmosphere / International Standard Atmosphere)
coesaRu = 8.31432 # J/mol/K

# Molecular mass of air as in COESA
# (US Atmosphere / International Standard Atmosphere)
coesaMair = 0.0289644 # kg/mol (mean molecular mass of air)