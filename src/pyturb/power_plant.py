from pyturb.control_volume import ControlVolume

class Intake(ControlVolume):
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
        print('solve intake')
        self.solved = True
        print(self.fluid.cp())
        return
        

class Combustor(ControlVolume):
    """
    """
    def __init__(self, stage):
        """
        """
        self.stage = stage
        return
        
        
    def solve(self):
        print('solve combustor')
        return