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
        products = ''
        hydrocarbon = False

        for element in self.fuel.thermo_prop.chemical_formula:
            if element is 'C':
                alpha = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "C{0:1.0f}".format(alpha) if not alpha==0 else self.stoichiometric_reaction
                hydrocarbon = True
            
            elif element is 'H':
                beta = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "H{0:1.0f}".format(beta) if not beta==0 else self.stoichiometric_reaction
                hydrocarbon = True if hydrocarbon else False
            
            elif element is 'O':
                gamma = self.fuel.thermo_prop.chemical_formula[element]
                reactants += "O{0:1.0f}".format(gamma) if not gamma==0 else self.stoichiometric_reaction

        oxidizer_fuel_ratio = alpha + beta/4 - gamma/2

        if self.oxidizer.gas_species is 'Air':
            delta = 1
            reactants += " + {0:1.3f}(O2 + 79/21 N2)".format(oxidizer_fuel_ratio)
            if hydrocarbon:
                products = "{0:1.0f} CO2 + {1:1.0f} H2O + {2:1.0f}·79/21 N2".format(alpha, beta/2, oxidizer_fuel_ratio)
            

        elif self.oxidizer.gas_species is 'O2':
            delta = 0
            reactants += " + {0:1.3f} O2".format(oxidizer_fuel_ratio)
            
            if hydrocarbon:
                products = "{0:1.0f} CO2 + {1:1.0f} H2O".format(alpha, beta/2)

        elif self.oxidizer.gas_species is 'O3':
            # XXX terminar
            hydrocarbon = True if hydrocarbon else False
        
        else:
            hydrocarbon = False
            raise NotImplementedError

        # Store stoichiometric coefficients:
        self._alpha = alpha
        self._beta = beta
        self._gamma = gamma
        self._delta = delta
        self._oxidizer_fuel_ratio = oxidizer_fuel_ratio
        

        if hydrocarbon:
            self.stoichiometric_reaction = reactants + ' --> ' + products


    @property
    def alpha(self):
        """
        """
        return self._alpha

        @property
    def beta(self):
        """
        """
        return self._beta

        @property
    def gamma(self):
        """
        """
        return self._gamma

        @property
    def delta(self):
        """
        """
        return self._delta

        @property
    def oxidizer_fuel_ratio(self):
        """
        """
        return self._oxidizer_fuel_ratio