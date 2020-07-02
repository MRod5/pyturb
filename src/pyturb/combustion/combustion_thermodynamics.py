"""
"""
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
import numpy as np
import warnings

oxidizers = ['Air', 'O2', 'O3', 'O2(L)', 'O3(L)']
fuels = ['hydrocarbon', 'C8H18,isooctane', '']

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
    def hcomb(self):
        """
        Combustion Enthalpy. [J/mol]
        """
        return self._hcomb


    def _classify_reactants(self):
        """
        """

        if not self.oxidizer.gas_species in self.oxidizer_list:
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

        alpha = 0
        beta = 0
        gamma = 0
        reactants = ""
        productsC = ""
        productsH = ""
        
        reactants_dict = {}
        products_dict = {}


        # TODO: Fuel mixtures will need a rework of alpha, beta, gamma coefficients

        # Reactants:
        for element in self.fuel.thermo_prop.chemical_formula:
            if element is "C":
                alpha = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "C{0:1.0f}".format(alpha) if not alpha==0 else self.stoichiometric_reaction

                # XXX: CO2 is formed if O2 is the oxidizer, not just when carbon appears
                productsC += "{0:1.0f}CO2".format(alpha)
            
            elif element is "H":
                beta = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "H{0:1.0f}".format(beta) if not beta==0 else self.stoichiometric_reaction

                # XXX: same. H2O is formed if oxygen is present.
                productsH += "{0:1.0f}H2O".format(beta/2)
            
            elif element is "O":
                gamma = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "O{0:1.0f}".format(gamma) if not gamma==0 else self.stoichiometric_reaction

        # Oxidizer to fuel ratio depending on C, H, O atoms
        self._oxidizer_fuel_ratio = alpha + beta/4 - gamma/2

        # Products:
        if not (productsC is "" and productsH is ""):
            products = productsC + " + " + productsH

        elif productsC is "":
            products = productsH

        else:
            products = productsC

        # If mixture has inert species
        if self.oxidizer.gas_species is "Air":
            delta = 1
            reactants += " + {0:1.3f}(O2 + 79/21 N2)".format(self.oxidizer_fuel_ratio)
            products += " + {0:1.3f}(79/21 N2)".format(self.oxidizer_fuel_ratio)

        else:
            delta = 0
            reactants += " + {0:1.3f} O2".format(self.oxidizer_fuel_ratio)

        # Combustion reaction:
        self._stoichiometric_reaction = reactants + " --> " + products
        self._reactants = reactants
        self._products = products
        
        # Coefficients:
        self._alpha = alpha
        self._beta = beta
        self._gamma = gamma
        self._delta = delta
        self._oxidizer_fuel_ratio = self.oxidizer_fuel_ratio

        # TODO: Must be an independent function, otherwise the FAR of a full reaction wont make sense
        # Fuel/Air ratio:
        if delta == 1:
            self._stoich_far = self.fuel.thermo_prop.Mg / (self.oxidizer_fuel_ratio/0.21*self.oxidizer.thermo_prop.Mg)
        
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

