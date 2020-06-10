from pyturb.fluid import Fluid

class Intake(Fluid):
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
        

class Combustor():
    """
    """
    def __init__(self, stage):
        """
        """
        self.stage = stage
        return
        
    def solve(self):
#            print('intake?', super().intake)
        print('solve combustor')
        return