"""

"""

import os
import numpy as np

nasa_thermoprop = r''
thermoprop_dir = os.path.dirname(__file__) 
nasa_thermoprop = os.path.join(thermoprop_dir, r'./nasa9.dat')

species='O2'

with open(nasa_thermoprop, 'r') as nasa_file:
    keep_searching = True
    
    while keep_searching:
        line = nasa_file.readline()
        print(line, line='')
        if line.split()[0]==species:
            # Get reference for the species data:
            data = nasa_file.readline()
            temp_intervals = int(data[0:2])
            Mg = float(data[52:65]) # g/mol
            deltaHfref = float(data[65:80])

            Temp_range = np.zeros([temp_intervals, 2])
            coeffs = np.zeros([temp_intervals, 7])
            for ii in range(temp_intervals):
                temprange_data = nasa_file.readline()
                coeffs_line_1 = nasa_file.readline()
                coeffs_line_2 = nasa_file.readline()

                Temp_range[ii][:]=[1,2]
                
                Temp_range[ii,:] = [float(temprange_data[0:11]), float(temprange_data[11:22])]
               
                coeffs[ii,:] = [float(coeffs_line_1[0:16].replace('D','E')), float(coeffs_line_1[16:32].replace('D','E')),
                                float(coeffs_line_1[32:48].replace('D','E')), float(coeffs_line_1[48:64].replace('D','E')),
                                float(coeffs_line_1[64:80].replace('D','E')), float(coeffs_line_2[0:16].replace('D','E')),
                                float(coeffs_line_2[16:32].replace('D','E'))]
                
            deltaHf0 = float(temprange_data[65:79])
            keep_searching = False
        elif line==None:
            keep_searching = False
        else:
            print(line)

            
print(keep_searching)
def cp(T):
    
    #cp_factor = coeffs
    #cp_gas = cp_factor * Rg
    return cp_gas
