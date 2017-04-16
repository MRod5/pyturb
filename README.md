# pyTurb
---

A python Gas Turbine package!

[TOC]

---

### Deploying code 
To deploy code run **setup.py**:
> python setup.py install

At the moment, only *numpy* and *scipy* are required to run.

### Contents
* **ISA.py:** International standard atmosphere [1]
* **air_model.py:** Thermally perfect or real gas approach for the air. For the real gas model, NASA and NACA models are proposed among others. [2], [3]
*  **isentropic_gas.py:** Isentropic gas relations, Mach number, stagnation enthalpy, and flow velocity.


### Examples
Examples are found at **./pyturb/notebooks**
1. isa_example.ipybn
2. air_model_example.ipynb
3. isentropic_gas_example.ipynb


---

### References
**[1]** - Getting to grips with aircraft performance. Airbus Operations, Customer Services.
**[2]** - *Equations tables and charts for compressible flow*. National Advisory Committee for Aeronautics, report 1135
**[3]** - *Coefficients for calculating thermodynamic and transport properties of individual species*. NASA Technical Memorandum 4513

April 2017