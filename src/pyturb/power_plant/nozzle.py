"""
Nozzle:
-------

Generic nozzle control volume. Extends from ControlVolume.

Implements 30 different thermodynamic properties and variables of the 
control volume.


MRodriguez 2020

"""

from pyturb.power_plant.control_volume import ControlVolume

class Nozzle(ControlVolume):
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
        print('solve nozzle')
        self.solved = True
        print(self.fluid.cp())
        return11