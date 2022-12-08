"""
combustion_thermodynamics:
--------------------------

Combustion reaction thermodynamical properties.

MRodriguez 2020

"""
from pyturb.gas_models import ThermoProperties
from pyturb.gas_models import PerfectIdealGas
from pyturb.gas_models import SemiperfectIdealGas
from pyturb.gas_models import GasMixture
from pyturb.utils import units
import numpy as np
import warnings

#TODO: Update gases:
oxidizers = ['Air', 'O', 'O2', 'O3'] # Allowed oxidizers

fuels = ['hydrocarbon',
         'CH4', 'C2H6', 'C3H8', 'C4H10', 'C5H12', 'C6H14', 'C7H16', 'C8H18',
         'C9H19', 'C10H8',
         'CH4O', 'CH3OCH3',
         'C2H2',
         'H2'] # Allowed fuels

inert_gases = ['He', 'Ar', 'N2',
               'CO2', 'CO',
               'H2O']

class Combustion(object):
    """
    Combustion:
    -----------

    This module is under development.

    Calculates the combustion between fuel and oxidizer.

    MRodriguez. 2020

    """

    def __init__(self, fuel, oxidizer):
        """
        Initiates a Combustion(). Fuel and Oxidizer must be provided.
        """

        if not(isinstance(fuel, PerfectIdealGas) or isinstance(fuel, SemiperfectIdealGas) or isinstance(fuel, IdealLiquid)):
            # Check the fluid is a Perfect or a Semiperfect gas fom pyturb
            raise TypeError("Object must be PerfectIdealGas, SemiperfectIdealGas or PerfectLiquid. Instead received {}".format(fuel))


        if not(isinstance(oxidizer, PerfectIdealGas) or isinstance(fuel, SemiperfectIdealGas)
                or isinstance(oxidizer, GasMixture) or isinstance(fuel, IdealLiquid)):
            # Check the fluid is a Perfect or a Semiperfect gas from pyturb
            raise TypeError("Object must be PerfectIdealGas, SemiperfectIdealGas or PerfectLiquid. Instead received {}".format(oxidizer))

        # Fuel and oxidizer lists:
        self.oxidizer_list = oxidizers
        self.fuel_list = fuels

        # Selected fuel and oxidizer for combustion reaction
        self.fuel = fuel
        self.oxidizer = oxidizer

        # Check if reactants are allowed fuels and oxidizers:
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
        Check fuel and oxidizer species. Oxidizer can be a gax mixtures with at least one of the allowed oxidizers.
        """

        # Check for the oxidizers:
        reactants_check = False
        if not (self.oxidizer.gas_species in self.oxidizer_list or self.oxidizer.gas_species == "mixture"):
            for oxidizer in self.oxidizer_list:
                if oxidizer in self.oxidizer.gas_species:
                    reactants_check = True
                    break

            if not reactants_check:
                warnings.warn("Requested oxidizer ({0}) not available. Available oxidizers: {1}".format(self.oxidizer.gas_species, self.oxidizer_list))
                return False

        # Check the mixture of oxidizers
        elif self.oxidizer.gas_species=="mixture":
            oxid_mix = self.oxidizer.mixture_gases['gas_species']
            for oxidizer in oxid_mix:
                if oxidizer in self.oxidizer_list:
                    reactants_check = True
                    break
            
            if not reactants_check:
                warnings.warn("Requested gas mixture ({0}) not available. Available oxidizers: {1}".format(oxid_mix, self.oxidizer_list))
                return False
    
        # Check fuels
        if not (self.fuel.gas_species in self.fuel_list or self.fuel.gas_species == "mixture"):
            for fuel in self.fuel_list:
                if fuel in self.fuel.gas_species:
                    reactants_check = True
                    break

            if not reactants_check:
                warnings.warn("Requested fuel ({0}) not available. Available fuels: {1}".format(self.fuel.gas_species, self.fuel_list))
                return False

        # Return with True, all reactants OK.
        return True


    def _combustion_stoichiometry_simple_reaction(self):
        """
        Stoichiometric reaction of a combustion with one mole of fuel and an oxidizer.
        """

        has_carbon = False
        has_hydrogen = False
        has_oxygen = False

        reactants = ""
        products = ""
        inerts = ""
        products_dictionary = {}
        reactants_dictionary = {}

        alpha = 0 # Carbons
        beta = 0  # Hydrogens
        gamma = 0 # Oxygens
        delta = 0 # Inerts//Nytrogen
       

        # TODO: Fuel mixtures will need a rework of alpha, beta, gamma coefficients

        ## Reactants:
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

        # Oxidizer
        if self.oxidizer.gas_species == "Air":
            # Air case, 21%O2 + 79%N2 (inert)
            delta = 1 #Has inert gas
            has_oxygen = True
            reactants += " + {0:1.3f}(O2 + 79/21 N2)".format(self.oxidizer_fuel_ratio)
            reactants_dictionary["Air"] = self.oxidizer_fuel_ratio

        elif self.oxidizer.gas_species == "O2" or self.oxidizer.gas_species == "O2(L)":
            delta = 0 # Oxicombustion

            for element in self.oxidizer.thermo_prop.chemical_formula:
                if element == "O":
                    has_oxygen = True
                    reactants += " + {0:1.3f} O{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])
                    reactants_dictionary["Air"] = self.oxidizer_fuel_ratio

        elif self.oxidizer.gas_species == "O3" or self.oxidizer.gas_species == "O3(L)":
            delta = 0 # Oxicombustion with ozone
            self._oxidizer_fuel_ratio = self.oxidizer_fuel_ratio * 2/3 # O3 instead of O2
            if element == "O":
                has_oxygen = True
                reactants += " + {0:1.3f} O{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])
            reactants_dictionary["O3"] = self.oxidizer_fuel_ratio

        elif self.oxidizer.gas_species == "O":
            delta = 0
            self._oxidizer_fuel_ratio = self.oxidizer_fuel_ratio * 2 # Monoatomic oxygen
            if element == "O":
                has_oxygen = True
                reactants += " + {0:1.3f} O{1:1.0f}".format(self.oxidizer_fuel_ratio, self.oxidizer.thermo_prop.chemical_formula[element])
            reactants_dictionary["O"] = self.oxidizer_fuel_ratio

        else:
            if self.oxidizer.gas_species == "mixture":
                a = b = c = d = 0 # Oxygen coefficient
                delta = 0
                for ii, gases in enumerate(self.oxidizer.mixture_gases['gas_species']):
                    if gases == "O":
                        has_oxygen = True
                        a = self.oxidizer.mixture_gases.loc[ii]['Ng'] # Monoatomic oxygen

                    elif gases == "O2":
                        has_oxygen = True
                        b = self.oxidizer.mixture_gases.loc[ii]['Ng'] # Diatomic oxygen

                    elif gases == "O3":
                        has_oxygen = True
                        c = self.oxidizer.mixture_gases.loc[ii]['Ng'] # Ozone


                self._oxidizer_fuel_ratio = self.oxidizer_fuel_ratio * 2 / (a + b*2 + c*3) # OFR is calculated per diatomic oxygen molecule
                
                if a!=0:
                    # O in mixture
                    reactants += " + {0:1.3f} O".format(a*self.oxidizer_fuel_ratio)
                    reactants_dictionary["O"] = a*self.oxidizer_fuel_ratio

                if b!=0:
                    # O2 in mixture
                    reactants += " + {0:1.3f} O2".format(b*self.oxidizer_fuel_ratio)
                    reactants_dictionary["O2"] = b*self.oxidizer_fuel_ratio

                if c!=0:
                    # O3 in mixture
                    reactants += " + {0:1.3f} O3".format(c*self.oxidizer_fuel_ratio)
                    reactants_dictionary["O3"] = c*self.oxidizer_fuel_ratio

                # Non-oxidizer gases in the mixture
                for ii, gases in enumerate(self.oxidizer.mixture_gases['gas_species']):
                    if gases in inert_gases:
                        delta = 1
                        d = self.oxidizer.mixture_gases.loc[ii]['Ng'] * self._oxidizer_fuel_ratio
                        inerts += " + {0:1.3f}{1}".format(d, gases)
                        reactants += inerts
                        reactants_dictionary[gases] = d


        ## Products:
        if has_oxygen:
            if has_carbon and not has_hydrogen: 
                # Carbon based fuel (non-hydrocarbon))
                products += "{0:1.3f}CO2".format(alpha)
                products_dictionary['CO2'] = alpha

            elif not has_carbon and has_hydrogen:
                # Hydrogen:
                products += "{0:1.3f}H2O".format(beta/2)
                products_dictionary['H2O'] = beta/2
            
            elif has_carbon and has_hydrogen:
                # Hydrocarbon case:
                products += "{0:1.3f}CO2 + {1:1.3f}H2O".format(alpha, beta/2)
                products_dictionary['CO2'] = alpha
                products_dictionary['H2O'] = beta/2

        #else:        
            # TODO: Complete with other fuels and oxidizers in the future

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


        # Stoichiometric Fuel/Air ratio or Fuel/Oxygen ratio:
        if self.oxidizer.gas_species == "Air":
            self._stoich_far = self.fuel.thermo_prop.Mg / (self.oxidizer_fuel_ratio/0.21*self.oxidizer.Mg)
            self._stoich_for = np.nan
        else:
            # TODO: FAR for non-air oxidizers
            self._stoich_far = np.nan
            self._stoich_for = 0
            for oxidizer in ['O', 'O2', 'O3']:
                iOxidizer = self.oxidizer.mixture_gases[self.oxidizer.mixture_gases['gas_species']==oxidizer]
                if not iOxidizer.empty:
                    self._stoich_for += (iOxidizer['Mg']*iOxidizer['Ng'])*reactants_dictionary[oxidizer]
            
            self._stoich_for = (reactants_dictionary[self.fuel.gas_species]*self.fuel.Mg) / self._stoich_for

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


    def stoich_adiabatic_flame_temp(self, T0):
        """
        Adiabatic flame temperatura for stoichiometric reaction:

        The adiabatic flame temperature is calculated assuming all the heat of the combustion is invested in
        heating the products and only the products, with no heat losses (adiabatic volume).

        Is the reaction is in stoichiometric FAR, then the adiabatic flame temperature is the highest temperature
        that can be obtained with the provided fuel and oxidizers.

        Calculated temperature is obtained in Kelvin [K]
        """
        T0_ = units.celsius_to_kelvin(T0)
        Tad = T0_ 
        
        self._Tad_st = Tad

        # TODO: Declare temperature of the reactants and of the combustor

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

    
    @property
    def Tad_st(self):
        """
        Adiabatic flame temperature for the stoichiometric reaction. [K]
        """
        return self._Tad_st


    @property
    def Tad(self):
        """
        Adiabatic flame temperature. [K]
        """
        return self._Tad