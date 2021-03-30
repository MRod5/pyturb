# pyTurb: A python Gas Turbine package!


**pyTurb** is a Python module for Thermodynamics, Combustion and Power Plants. In **pyTurb** you will find:
- ISA model (Standard Atmosphere)
- Thermodynamic properties of gases (NASA Glenn Coefficients)
- Perfect and Semiperfect Ideal Gas relations and Gas Mixtures
- Isentropic Flow relations
- Combustion properties and kinetics
- Power Plant control volumes (nozzle, intake, combustor...)
- Propulsion equations (jet engine performance)

You will wind examples about how to use pyTurb in the *Notebooks* folder.

pyTurb is currently under development. In the near future, compressor, turbine and more combustion features will be added, as well as jet engine performance equations.


---
## How to use pyTurb:


The easiest way to use **pyTurb** is to install it with pip and the `setup.py` file. Download **pyTurb** and then, from the source folder, type in a command window:

>`pip install .`

If you want to have **pyTurb** in editable mode:

>`pip install -e .`

With pyTurb installed, you can import any *Gas Model*, *Combustion thermodynamics* or *Power Plant* features. You can find use-cases and examples about how to execute pyTurb at the `notebooks` folders.

## Notebooks with pyTurb examples:

You will find different Notebooks with examples about to use **pyTurb**:
- `ISA examples.ipynb`: Standard atmosphere example (temperature, pressure, density and altitude for different layers).
- `Perfect and Semiperfect gas models.ipynb`: gas constant, specific heat capacity, heat ratio etc. Examples as perfect gas (constant values) or semiperfect gas (as a function of the temperature).
- `Gas Mixtures`: Mixtures of different pure substances. May be a perfect (constant values) or semiperfect (as a function of temperature) gas.
- `Isentropic flow relations.ipynb`: Examples about basic isentropic flow relations: stagnation temperature, stagnation pressue, Mach number, sound speed, free stream velocity...
- ...


## Control Volumes:

`src\pyturb\power_plant\`

- Intake (`intake.py`): Diffuser control volume (lowers kinetic energy, raises sensible enthalpy)
- Nozzle (`nozzle.py`): Nozzle control volume (raises kinetic energy, lowers, T and p)
- Combustor (`combustor.py`): Heat power provided a fuel and oxidizer
- more to come...


## Gas Models:

`src\pyturb\gas_models\`

- *Perfect Ideal Gas* and *Semiperfect Ideal Gas* approaches (`perfect_ideal_gas.py`, `semiperfect_ideal_gas.py`):
  - Perfect gas: Constant cp, cv, gamma
  - Semiperfect gas: Specific heats are functions of the temperature. cv(T), cp(T), gamma(T)
  - Ideal gas: pv=RgT, Rg=cp-cv
- *Gas Mixtures* of ideal gases (`gas_mixture.py`):
  - Mixture properties (molar fraction, mass fraction, cp, cv, gamma...)
  - Mixture may be *Perfect* or *Semiperfect*
- Combustion thermodynamics (`combustion_thermodynamics.py`):
  - Stoichiometric coefficients, LHV, HHV, heat of combustion
  - Fuel mixtures and oxidizer mistures
  - Fuel-air ratio
- Implemented the International Standard Atmosphere (COESA 1975) from 0m to 84852m (`isa.py`)


## Utils:

`src\pyturb\utils\`

- Conversion of units to and from international system (`units.py`)
- Universal constants (`constants.py`)
- Numerical iterators: root finder with variable step and equation iterator for 2 independent variables (`numerical_iterators.py`)


---

## References

**[1]** - *Equations tables and charts for compressible flow*. National Advisory Committee for Aeronautics, report 1135

**[2]** - *Coefficients for calculating thermodynamic and transport properties of individual species*. NASA Technical Memorandum 4513

**[3]** - *NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species*. NASA / TP - 2002-211556

**[4]** - *Defining constants, equations, and abbreviated tables of the 1975 U.S. Standard Atmosphere*. NASA Technical Report TR R-459

---

Marcos Rodríguez

2020.
