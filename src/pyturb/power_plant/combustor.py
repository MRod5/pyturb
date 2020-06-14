"""
Combustor:
----------

Generic combustion chamber control volume. Extends from ControlVolume.

Implements 30 different thermodynamic properties and variables of the 
control volume.


MRodriguez 2020

"""

from pyturb.power_plant.control_volume import ControlVolume

class Combustor(ControlVolume):
    """
    """
    def __init__(self, stage, fluid):
        """
        """
        self.stage = stage
        super().__init__(fluid)
            
        return
        
    def solve(self):
        """
        """
        print('solve cc')
        self.solved = True
        print(self.fluid.cp())
        return