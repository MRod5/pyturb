"""
Intake:
-------

Generic intake (diffuser) control volume. Extends from ControlVolume.

Implements 30 different thermodynamic properties and variables of the 
control volume.


MRodriguez 2020

"""

from pyturb.power_plant.control_volume import ControlVolume

class Intake(ControlVolume):
    """
    """
    def __init__(self, stage, fluid):
        """
        """
        self.stage = stage
        super().__init__(fluid)
            
        return None
        
    
    def inputs(self, pe, Te, ve, adiab_efficiency=1):
        """
        Set basic inputs of a generic intake (diffuser):
            + Static pressure, static temperature and flow velocity at the entrance
            + Adiabatic efficiency of the diffuser
        """

        self.p_e = pe
        self.T_e = Te
        self.vel_e = ve
        
        return 

    def solve(self):
        """
        """
        print('solve intake')
        self.solved = True
        print(self.fluid.cp())
        return