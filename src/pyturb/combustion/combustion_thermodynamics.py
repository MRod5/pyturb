"""
combustion_thermodynamics:
--------------------------

Combustion reaction thermodynamical properties.

MRodriguez 2020

"""
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
import numpy as np
import warnings

oxidizers = ['Air', 'O', 'O2', 'O3', 'O2(L)', 'O3(L)']

fuels = ['hydrocarbon', 'C8H18,isooctane', 
         'CH4', 'C2H6', 'C3H8', 'C4H10', 'C5H12', 'C6H14', 'C7H16', 'C8H18',
         'CH4O', 'C2H6O']

inert_gases = ['He', 'Ar', 'N2',
               'CO2', 'CO']

class Combustion(object):
    """
    """

    def __init__(self, fuel, oxidizer):
        """
        """

        if not(isinstance(fuel, PerfectIdealGas) or isinstance(fuel, SemiperfectIdealGas) or isinstance(fuel, IdealLiquid)):
            # Check the flow is a Perfect or a Semiperfect gas fom pyturb
            raise TypeError("Object must be PerfectIdealGas, SemiperfectIdealGas or PerfectLiquid. Instead received {}".format(fluid))


        if not(isinstance(fuel, PerfectIdealGas) or isinstance(fuel, SemiperfectIdealGas) or isinstance(fuel, IdealLiquid)):
            # Check the flow is a Perfect or a Semiperfect gas from pyturb
            raise TypeError("Object must be PerfectIdealGas, SemiperfectIdealGas or PerfectLiquid. Instead received {}".format(fluid))

        
        self.oxidizer_list = oxidizers
        self.fuel_list = fuels


        self.fuel = fuel
        self.oxidizer = oxidizer


        reactants_status = self._classify_reactants()

        if not reactants_status:
            raise ValueError("Unknown fuel and oxidizer")
        else:
            self._alpha = 0
            self._beta = 0
            self._gamma = 0
            self._delta = 0


        return


    # Class properties for combustion thermodynamics
    @property
    def reactants(self):
        """
        Reactants of the combustion reaction.
        """
        return self._reactants


    @property
    def products(self):
        """
        Products of the combustion reaction.
        """
        return self._products


    @property
    def stoichiometric_reaction(self):
        """
        Stoichiometric reaction of the combustion reaction.
        """
        return self._stoichiometric_reaction


    @property
    def alpha(self):
        """
        Moles of carbon present in the fuel molecule, per mole of fuel.
        """
        return self._alpha

    
    @property
    def beta(self):
        """
        Moles of hydrogen present in the fuel molecule, per mole of fuel.
        """
        return self._beta


    @property
    def gamma(self):
        """
        Moles of oxygen present in the fuel molecule, per mole of fuel.
        """
        return self._gamma


    @property
    def delta(self):
        """
        Nytrogen present in the oxyder. 1 if true, 0 if false.
        """
        return self._delta


    @property
    def oxidizer_fuel_ratio(self):
        """
        Oxidizer to fuel stoichiometric molar ratio.
        """
        return self._oxidizer_fuel_ratio


    @property
    def stoich_far(self):
        """
        Stoichiometric fuel-air ratio.
        """
        return self._stoich_far

    
    @property
    def LHV(self):
        """
        Lower heating value. [J/kg]
        """
        return self._LHV


    @property
    def HHV(self):
        """
        Higher heating value. [J/kg]
        """
        return self._HHV


    @property
    def hcomb_l(self):
        """
        Combustion Enthalpy, considering all products are condensed. [J/mol]
        """
        return self._hcomb_l


    @property
    def hcomb_g(self):
        """
        Combustion Enthalpy, considered all products are vaporized. [J/mol]
        """
        return self._hcomb_g


    def _classify_reactants(self):
        """
        Check fuel and oxidizer species.
        """

        if not (self.oxidizer.gas_species in self.oxidizer_list or self.oxidizer.gas_species == "mixture"):
            warnings.warn("Requested oxidizer ({0}) not available. Available oxidizers: {1}".format(self.oxidizer.species, self.oxidizer_list))
            return False
    
        if not self.fuel.gas_species in self.fuel_list:
            warnings.warn("Requested fuel ({0}) not available. Available fuels: {1}".format(self.fuel.species, self.fuel_list))
            return False

        return True


    def _combustion_stoichiometry_simple_reaction(self):
        """
        Stoichiometric reaction of a combustion with one molecule of fuel and an oxidizer.
        """

        has_carbon = False
        has_hydrogen = False
        has_oxygen = False

        reactants = ""
        products = ""
        inerts = ""
        products_dictionary = {}
        reactants_dictionary = {}

        alpha = 0
        beta = 0
        gamma = 0
        delta = 0
       

        # TODO: Fuel mixtures will need a rework of alpha, beta, gamma coefficients

        # Reactants:
        reactants_dictionary[self.fuel.gas_species] = 1
        for element in self.fuel.thermo_prop.chemical_formula:
            if element == "C":
                alpha = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "C{0:1.0f}".format(alpha) if not alpha==0 else reactants
                has_carbon = True

            elif element == "H":
                beta = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "H{0:1.0f}".format(beta) if not beta==0 else reactants
                has_hydrogen = True

            
            elif element == "O":
                gamma = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "O{0:1.0f}".format(gamma) if not gamma==0 else reactants
                has_oxygen = True
        

        # Oxidizer to fuel ratio depending on C, H, O atoms, assuming O2 as oxidizer
        self._oxidizer_fuel_ratio = alpha + beta/4 - gamma/2

        # TODO: Rework products. Mixtures must be accepted...
        # Oxidizer
        if self.oxidizer.gas_species == "Air":
            delta = 1
            has_oxygen = True
            reactants += " + {0:1.3f}(O2 + 79/21 N2)".format(self.oxidizer_fuel_ratio)

        elif self.oxidizer.gas_species == "O2" or self.oxidizer.gas_species == "O2(L)":
            delta = 0
            for element in self.oxidizer.thermo_prop.chemical_formula:
                if element == "O":
                    has_oxygen = True
                    reactants += " + {0:1.3f} O{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])

                if element == "N":
                    delta = 1
                    reactants += " + {0:1.3f} N{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])

        elif self.oxidizer.gas_species == "O3" or self.oxidizer.gas_species == "O3(L)":
            self._oxidizer_fuel_ratio = self.oxidizer_fuel_ratio * 2/3
            if element == "O":
                    has_oxygen = True
                    reactants += " + {0:1.3f} O{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])

            if element == "N":
                delta = 1
                reactants += " + {0:1.3f} N{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])

        elif self.oxidizer.gas_species == "O":
            self._oxidizer_fuel_ratio = self.oxidizer_fuel_ratio * 2
            if element == "O":
                has_oxygen = True
                reactants += " + {0:1.3f} O{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])

            if element == "N":
                delta = 1
                reactants += " + {0:1.3f} N{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])


        else:
            if self.oxidizer.gas_species == "mixture":
                a = b = c = d = 0
                delta = 0
                for ii, gases in enumerate(self.oxidizer.mixture_gases['gas_species']):
                    if gases == "O":
                        has_oxygen = True
                        a = self.oxidizer.mixture_gases.loc[ii]['Ng']

                    elif gases == "O2":
                        has_oxygen = True
                        b = self.oxidizer.mixture_gases.loc[ii]['Ng']

                    elif gases == "O3":
                        has_oxygen = True
                        c = self.oxidizer.mixture_gases.loc[ii]['Ng']


                self._oxidizer_fuel_ratio = self.oxidizer_fuel_ratio * 2 / (a + 2*b + 3*c)
                if a!=0:
                    reactants += " + {0:1.3f} O".format(a*self.oxidizer_fuel_ratio)
                if b!=0:
                    reactants += " + {0:1.3f} O2".format(b*self.oxidizer_fuel_ratio)
                if c!=0:
                    reactants += " + {0:1.3f} O3".format(c*self.oxidizer_fuel_ratio)

                for ii, gases in enumerate(self.oxidizer.mixture_gases['gas_species']):
                    if gases in inert_gases:
                        delta = 1
                        d = self.oxidizer.mixture_gases.loc[ii]['Ng'] * self._oxidizer_fuel_ratio
                        inerts += " + {0:1.0f}{1}".format(d, gases)
                        reactants += inerts


        # Products:
        if has_oxygen:
            # Hydrocarbon and hydorgen case:
            if has_carbon and not has_hydrogen: 
                products += "{0:1.0f}CO2".format(alpha)
                
                if "CO2" in products_dictionary.keys():
                    products_dictionary['CO2'] = products_dictionary['CO2'] + alpha
                else:
                    products_dictionary['CO2'] = alpha

            elif not has_carbon and has_hydrogen:
                products += "{0:1.0f}H2O".format(beta/2)

                if "H2O" in products_dictionary.keys():
                    products_dictionary['H2O'] = products_dictionary['CO2'] + beta/2
                else:
                    products_dictionary['H2O'] = beta/2
            
            elif has_carbon and has_hydrogen:
                products += "{0:1.0f}CO2 + {1:1.0f}H2O".format(alpha, beta/2)
                if "CO2" in products_dictionary.keys():
                    products_dictionary['CO2'] = products_dictionary['CO2'] + alpha
                else:
                    products_dictionary['CO2'] = alpha
                
                if "H2O" in products_dictionary.keys():
                    products_dictionary['H2O'] = products_dictionary['CO2'] + beta/2
                else:
                    products_dictionary['H2O'] = beta/2

        #else:        
            # TODO: Complete with other ozidizers in the future

        # Nytrogen present:        
        if delta == 1:
            if self.oxidizer.gas_species == "Air":
                products += " + {0:1.3f}(79/21 N2)".format(self.oxidizer_fuel_ratio)
            else:
                products += inerts

       
        # Combustion reaction:
        self._stoichiometric_reaction = reactants + " --> " + products
        self._reactants = reactants
        self._products = products
        

        # Coefficients:
        self._alpha += alpha
        self._beta += beta
        self._gamma += gamma
        self._delta = delta if delta == 1 else self._delta


        # TODO: Must be an independent function, otherwise the FAR of a full reaction wont make sense
        # Fuel/Air ratio:
        if self.oxidizer.gas_species == "Air":
            self._stoich_far = self.fuel.thermo_prop.Mg / (self.oxidizer_fuel_ratio/0.21*self.oxidizer.Mg)
        
        else:
            self._stoich_far = np.nan


    def combustion_stoichiometry(self):
        """
        Stoichiometric reaction of the combustion process. If the fuel is a mixture,
        the combustion is decomposed in simple reactions with one molecule of fuel
        each. The global reaction stoichiometry is calculated adding up the 
        stoichiometric coefficients of the simpler reaction.
        """
        # TODO: This function should classify the fuel. If the fuel is a simple molecule,
        # TODO: call _combustion_stoichiometry_simple_reaction, otherwise identify all
        # TODO: molecules in the fuel mixture and then call the simple reaction
        self._combustion_stoichiometry_simple_reaction()
        return
        
    
    def heat_of_combustion(self):
        """
        """

        # TODO: calculate heat of combustion as -Qp (heat of comburtion by definition). If water is formed, 
        # calculate both condensed and vapor water

        # TODO: With the heat of combustion calculated, obtain the HHV and LHV

        