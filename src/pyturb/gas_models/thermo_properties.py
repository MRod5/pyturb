"""
Thermodynamic properties:
-------------------------

Thermodynamic properties of different species (NASA Glenn coefficients).
Creates an object with
    
    
MRodriguez. 2020

"""

import os
import numpy as np

class ThermoProperties():
    def __init__(self, species=None):
        """
        """
        
        thermoprop_dir = os.path.dirname(__file__) 
        self.nasa_thermoprop = os.path.join(thermoprop_dir, r'./nasa9.dat')
        self.species_list = self.get_species_name()

        if species is not None:
            if type(species) is str:
                self.species = species
            
                self.get_species_data()
            else:
                raise TypeError("species must be type str, instead {} was provided".format(type(species)))
        else:
            return None
        
        return None
        
        
    def __str__(self):
        printable_thermoprop = "Species: {0}\tMg={1} g/mol\tdeltaHf_ref={2} J/mol\tdeltaHf_0K={3} J/mol\n".format(self.species, self.Mg, self.deltaHf_ref, self.deltaHf_0K)
        printable_thermoprop += "-->Chemical formula: {}\n".format(self.chemical_formula)
        for ii in range(self.temp_intervals):
            printable_thermoprop += "-->Tinterval: [{0}:{1}] K\n".format(self.temp_range[ii,0], self.temp_range[ii,1])
            printable_thermoprop += "   Coefficients:  order -2  |  order -1  |  order 0  |  order 1  |  order 2  |  order 3  |  order 4\n".format()
            printable_thermoprop += "                 {0:10.3e}   {1:10.3e}  {2:10.3e}  {3:10.3e}  {4:10.3e}  {5:10.3e}  {6:10.3e}\n".format(*self.coefficients[ii])
            printable_thermoprop += "\n"
        
        return printable_thermoprop


    def get_species_data(self):
        with open(self.nasa_thermoprop, 'r') as nasa_file:
            keep_searching = True
            self.Mg = None
            self.deltaHf_ref = None
            self.temp_range = None
            self.temp_intervals = None
            self.coefficients = None
            self.deltaHf_0K = None
            self.chemical_formula = {}
            

            while keep_searching:
                line = nasa_file.readline()
                
                if line=='':
                    keep_searching = False
                    
                    if self.Mg == None:
                        raise ValueError("Unknown species: {}".format(self.species))


                elif line.split()[0]==self.species:
                    # Get reference for the species data:
                    data = nasa_file.readline()
                    self.temp_intervals = int(data[0:2])
                    self.Mg = float(data[52:65])
                    self.deltaHf_ref = float(data[65:80])
                    formula = data[10:50]
                    
                    for i_element in range(5):
                        element = formula[8*(i_element):8*(i_element+1)]
                        atom = element[0:2]
                        if not atom=='  ':
                            n_atoms = float(element[2:])
                            atom = atom.strip()
                            self.chemical_formula[atom] = n_atoms

                    self.temp_range = np.zeros([self.temp_intervals, 2])
                    self.coefficients = np.zeros([self.temp_intervals, 7])
                    self.integration_cts = np.zeros([self.temp_intervals, 2])
                    
                    for ii in range(self.temp_intervals):
                        temprange_data = nasa_file.readline()
                        coeffs_line_1 = nasa_file.readline()
                        coeffs_line_2 = nasa_file.readline()

                        self.temp_range[ii][:]=[1,2]
                        self.temp_range[ii,:] = [float(temprange_data[0:11]), float(temprange_data[11:22])]
                        self.coefficients[ii,:] = [float(coeffs_line_1[0:16].replace('D','E')), float(coeffs_line_1[16:32].replace('D','E')),
                                        float(coeffs_line_1[32:48].replace('D','E')), float(coeffs_line_1[48:64].replace('D','E')),
                                        float(coeffs_line_1[64:80].replace('D','E')), float(coeffs_line_2[0:16].replace('D','E')),
                                        float(coeffs_line_2[16:32].replace('D','E'))]
                        
                        self.integration_cts[ii,:] = [float(coeffs_line_2[48:64].replace('D','E')), float(coeffs_line_2[64:80].replace('D','E'))]
                        

                    self.deltaHf_0K = float(temprange_data[65:79])
                    keep_searching = False
                    
                        
                else:
                    pass

                
    def get_species_name(self):
        """
        """
        self.species_list = []

        with open(self.nasa_thermoprop, 'r') as nasa_file:
            keep_searching = True

            while keep_searching:
                line = nasa_file.readline()

                if line=='':
                    keep_searching = False
                else:
                    try:
                        float(line.split()[0][:2])
                    except:
                        self.species_list.append(line.split()[0])

        return self.species_list

    
    def is_available(self, species):
        """
        """

        if type(species) is str:
            if species in self.species_list:
                return True
            else:
                return False
                
        else:
            raise TypeError('Gas species must be of type str, instead {} was provided'.format(type(species)))