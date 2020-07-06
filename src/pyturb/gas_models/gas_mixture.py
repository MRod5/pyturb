"""
"""
import pyturb.utils.constants as cts
from pyturb.gas_models.gas import Gas
from pyturb.gas_models.perfect_ideal_gas import PerfectIdealGas
from pyturb.gas_models.semiperfect_ideal_gas import SemiperfectIdealGas
import numpy as np
import pandas as pd
import warnings


class GasMixture(object):
    """
    """

    def __init__(self, gas_model = "perfect", mixture=None):
        """
        """

        # Gas model selector:
        if gas_model.lower() == "perfect":
            self.gas_model = PerfectIdealGas
        elif gas_model.lower() == "semi-perfect" or gas_model.lower() == "semiperfect":
            self.gas_model = SemiperfectIdealGas
        else:
            self.gas_model = None
            raise ValueError("gas_model may be 'perfect' or 'semi-perfect', instead received {}".format(gas_model))
        
        self._n_species = 0
        self._mixture_gases_columns = ['gas_species', 'gas_properties', 'Ng', 'Mg', 'mg', 'Rg', 'molar_frac', 'mass_frac']
        self._mixture_gases = pd.DataFrame(columns=self._mixture_gases_columns)


        # Mixture initialization:
        if not mixture is None:
            if not type(mixture) is dict:
                warnings.warn("mixture must be a dictionary with keys:=species and value:=moles. Instead received {}".format(mixture))
            else:
                #TODO: call add_gas for each pure substance in the mixture
                print("")

        return



    def add_gas(self, species, moles=None, mass=None):
        """
        """

        if moles is None and mass is None:
            raise ValueError('Quantity (moles or mass) of gas must be specified in add_gas.')
        else:

            pure_substance = self.gas_model(species)

            if mass is None:
                moles_ = moles # mol
                mass_ = moles_ * pure_substance.Rg * 1e-3 # kg

            elif moles is None:
                mass_ = mass #kg
                moles_ = mass_ / pure_substance.Rg * 1e3 # mol

            else:
                warnings.warn("mass ({0}kg) will be dismised and recalculated with the moles quantity provided: {1}mol.".format(mass, moles))
                moles_ = moles
                mass_ = mass_ = moles_ * pure_substance.Rg * 1e3 # kg



        subst_props = {'gas_species': pure_substance.gas_species, 'gas_properties': pure_substance, 'Ng': moles_, 'Mg': pure_substance.thermo_prop.Mg, 'mg': mass_, 'Rg': pure_substance.Rg, 'molar_frac': np.nan, 'mass_frac': np.nan}

        new_gas = len(self._mixture_gases.index)
        self._mixture_gases.loc[new_gas] = subst_props

        self._update_mixture_properties()

        return


    def _update_mixture_properties(self):
        """
        """

        self._n_species = len(self.mixture_gases.index)

        self._Ng = np.sum(self.mixture_gases['Ng'])
        self._mg = np.sum(self.mixture_gases['mg'])

        self._mixture_gases['molar_frac'] = self.mixture_gases['Ng'] / self.Ng
        self._mixture_gases['mass_frac'] = self.mixture_gases['mg'] / self.mg

        self._Mg = 0
        for ii, xi in enumerate(self.mixture_gases['molar_frac']):
            self._Mg += xi * self.mixture_gases.loc[ii]['Mg']


        return    


    @property
    def n_species(self):
        """
        """
        return self._n_species


    @property
    def mixture_gases(self):
        """
        """
        return self._mixture_gases
    
    
    @property
    def Ng(self):
        """
        """
        return self._Ng
    

    @property
    def mg(self):
        """
        """
        return self._mg
    

    @property
    def Ru(self):
        """
        Get the Ideal Gas Law constant Ru [J/mol/K]
        """
        
        Ru = cts.Ru
        return Ru

    
    @property
    def Rg(self):
        """
        Get the Mixture Gas constant Rg =  Ru/Mg [J/kg/K]
        """
        
        Rg = self.Ru/self.Mg*1e3
        return Rg


    @property
    def Mg(self):
        """
        Get the Mixture molecular mass [g/mol]
        """
        return self._Mg

