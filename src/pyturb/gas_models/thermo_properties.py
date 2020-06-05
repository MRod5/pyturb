"""

"""

import os
import numpy as np

class ThermoProperties():
    def __init__(self, species):
        thermoprop_dir = os.path.dirname(__file__) 
        self.nasa_thermoprop = os.path.join(thermoprop_dir, r'./nasa9.dat')
        self.species = species
        
        self.get_species_data()
        
        return


    def get_species_data(self):
        with open(self.nasa_thermoprop, 'r') as nasa_file:
            keep_searching = True
            self.Mg = None
            self.deltaHf_ref = None
            self.temp_range = None
            self.temp_intervals = None
            self.coefficients = None
            self.deltaHf_0K = None
            

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

                    self.temp_range = np.zeros([self.temp_intervals, 2])
                    self.coefficients = np.zeros([self.temp_intervals, 7])
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

                    self.deltaHf_0K = float(temprange_data[65:79])
                    keep_searching = False
                    
                        
                else:
                    pass

                
