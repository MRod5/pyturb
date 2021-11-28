"""
combustion_thermodynamics:
--------------------------

Combustion reaction thermodynamical properties.

MRodriguez 2020

"""
from pyturb.gas_models.thermo_properties import ThermoProperties
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
import numpy as np
import warnings

#TODO: Update gases:
oxidizers = ['Air', 'O', 'O2', 'O3', 'O2(L)', 'O3(L)']

fuels = ['hydrocarbon', 'C8H18,isooctane', 
         'CH4', 'C2H6', 'C3H8', 'C4H10', 'C5H12', 'C6H14', 'C7H16', 'C8H18',
         'CH4O', 'CH3OCH3',
         'H2']

inert_gases = ['He', 'Ar', 'N2',
               'CO2', 'CO']

class Combustion(object):
    """
    Combustion:
    -----------


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


    def _classify_reactants(self):
        """
        Check fuel and oxidizer species.
        """

        if not (self.oxidizer.gas_species in self.oxidizer_list or self.oxidizer.gas_species == "mixture"):
            warnings.warn("Requested oxidizer ({0}) not available. Available oxidizers: {1}".format(self.oxidizer.gas_species, self.oxidizer_list))
            return False
    
        if not (self.fuel.gas_species in self.fuel_list or self.fuel.gas_species == "mixture"):
            warnings.warn("Requested fuel ({0}) not available. Available fuels: {1}".format(self.fuel.gas_species, self.fuel_list))
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
        # Stoichiometric combustion are calculated per unit mole of fuel. Quantity is adjusted afterwards.
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
            reactants_dictionary["Air"] = self.oxidizer_fuel_ratio

        elif self.oxidizer.gas_species == "O2" or self.oxidizer.gas_species == "O2(L)":
            delta = 0
            for element in self.oxidizer.thermo_prop.chemical_formula:
                if element == "O":
                    has_oxygen = True
                    reactants += " + {0:1.3f} O{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])
                    reactants_dictionary["Air"] = self.oxidizer_fuel_ratio

        elif self.oxidizer.gas_species == "O3" or self.oxidizer.gas_species == "O3(L)":
            self._oxidizer_fuel_ratio = self.oxidizer_fuel_ratio * 2/3
            if element == "O":
                has_oxygen = True
                reactants += " + {0:1.3f} O{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])
            reactants_dictionary["O3"] = self.oxidizer_fuel_ratio

        elif self.oxidizer.gas_species == "O":
            self._oxidizer_fuel_ratio = self.oxidizer_fuel_ratio * 2
            if element == "O":
                has_oxygen = True
                reactants += " + {0:1.3f} O{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])
            reactants_dictionary["O"] = self.oxidizer_fuel_ratio

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
                    reactants_dictionary["O"] = a*self.oxidizer_fuel_ratio
                if b!=0:
                    reactants_dictionary["O2"] = b*self.oxidizer_fuel_ratio
                if c!=0:
                    reactants += " + {0:1.3f} O3".format(c*self.oxidizer_fuel_ratio)
                    reactants_dictionary["O3"] = c*self.oxidizer_fuel_ratio


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
                products_dictionary['CO2'] = alpha

            elif not has_carbon and has_hydrogen:
                products += "{0:1.0f}H2O".format(beta/2)
                products_dictionary['H2O'] = beta/2
            
            elif has_carbon and has_hydrogen:
                products += "{0:1.0f}CO2 + {1:1.0f}H2O".format(alpha, beta/2)
                products_dictionary['CO2'] = alpha
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
        self._reactants_dictionary_monofuel = reactants_dictionary
        self._products_dictionary_monofuel = products_dictionary
        

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

        if self.fuel.gas_species == "mixture":
            # TODO
            print("mixture")
        else:
            self._combustion_stoichiometry_simple_reaction()
            self._reactants_dictionary = self._reactants_dictionary_monofuel
            self._products_dictionary = self._products_dictionary_monofuel

        return
        
    
    def heat_of_combustion(self):
        """
        Calculates the heat of a combustion under the model of constant pressure
        combustor:
            -qp = Sum_i(N_i*h_fi) - Sum_j(N_j*h_fj) [J/mol]
        Heat of combustion is calculated per unit of mole of fuel and then it is 
        scaled to the mole quantity of fuel present in the reaction.

        If the products contain water, the heat of combustion is calculated assuming
        all water is condensed (l) and vaporized (g), providing two different values
        of the heat of combustion:
            - hcomb_g: heat of combustion with H2O(g) vaporized [J/mol]
            - hcomb_l: heat of combustion with H2O(l) condensed [J/mol]

        With the heat of combustion per unit mole, the Higher Heating Value and Lower
        Heating Value may be calculated:
            HHV = LHV + NH2O*MH2O*hfg / (Nfuel*Mfuel) [MJ/kg]
            LHV = hcomb_g / Mfuel * 1e3 g/kg * 1e-6 MJ/J [MJ/kg]
        """
        
        # Loop for calculating the formation enthalpies of the reactants:
        qp_r = 0
        reactants = self.reactants_dictionary.keys()
        
        for element in reactants:
            if element == "Air":
                # Air is excluded from the thermodynamic_properties call, its heat of formation is negligible
                qp_r += 0
            else:
                thermo_prop = ThermoProperties(element)
                qp_r += thermo_prop.deltaHf_ref * self.reactants_dictionary[element]
        

        # Loop for calculating the formation enthalpies of the products
        products = self.products_dictionary.keys()

        qp_pl = 0
        qp_pg = 0
        for element in products:
            print(element)
            if element == "Air":
                # Air is excluded from the thermodynamic_properties call, its heat of formation is negligible
                qp_pl += 0
                qp_pg += 0
            else:
                if element == 'H2O':
                    thermo_prop = ThermoProperties('H2O(L)')
                    qp_pl += thermo_prop.deltaHf_ref * self.products_dictionary[element]

                    thermo_prop = ThermoProperties('H2O')
                    qp_pg += thermo_prop.deltaHf_ref * self.products_dictionary[element]
                else:
                    thermo_prop = ThermoProperties(element)
                    qp_pl += thermo_prop.deltaHf_ref * self.products_dictionary[element]
                    qp_pg += thermo_prop.deltaHf_ref * self.products_dictionary[element]
        
        
        self._hcomb_g = qp_r - qp_pg
        self._hcomb_l = qp_r - qp_pl

        self._HHV = self.hcomb_l / self.fuel.Mg * 1e3
        self._LHV = self.hcomb_g / self.fuel.Mg * 1e3

        # TODO: calculate heat of combustion as -Qp (heat of combustion by definition). If water is formed, 
        # calculate both condensed and vapor water

        # TODO: With the heat of combustion calculated, obtain the HHV and LHV

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


    # Class properties for combustion thermodynamics
    @property
    def reactants_dictionary(self):
        """
        Reactants dictionary [gas_species]: moles
        """
        return self._reactants_dictionary


    @property
    def products_dictionary(self):
        """
        Products dictionary [gas_species]: molesProducts of the combustion reaction.
        """
        return self._products_dictionary


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