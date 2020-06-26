"""
"""
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
import warnings

oxidizers = ['Air', 'O2', 'O3', 'O2(L)', 'O3(L)']
fuels = ['hydrocarbon', 'C8H18,isooctane', '']

class ThermoCombustion(object):
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

    def stoichiometry(self):
        """
        """

        alpha = 0
        beta = 0
        gamma = 0
        self.stoichiometric_reaction = ''
        reactants = ''
        product = ''
        print(self.fuel.thermo_prop.chemical_formula)
        for element in self.fuel.thermo_prop.chemical_formula:
            if element is "C":
                alpha = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "C{0:1.0f}".format(alpha) if not alpha==0 else self.stoichiometric_reaction
            
            elif element is "H":
                beta = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "H{0:1.0f}".format(beta) if not beta==0 else self.stoichiometric_reaction
            
            elif element is "O":
                gamma = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "O{0:1.0f}".format(gamma) if not gamma==0 else self.stoichiometric_reaction

        oxidizer_fuel_ratio = alpha + beta/4 - gamma/2

        if self.oxidizer.gas_species is "Air":
            delta = 1

            reactants += " + {0:1.3f}(O2 + 79/21 N2)".format(oxidizer_fuel_ratio)

        else:
            delta = 0
            
            reactants += " + {0:1.3f} O2".format(oxidizer_fuel_ratio)

        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.delta = delta

