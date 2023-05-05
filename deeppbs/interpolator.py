# third party packages

class Interpolator(object):
    def __init__(self, fileName):
        try:
            from gridData import Grid
        except ModuleNotFoundError:
            raise ModuleNotFoundError("The dependency 'GridDataFormats' is required for this functionality!")
        
        self.grid = Grid(fileName) # stores the grid data
            
    def __call__(self, xyz):
        return self.grid.interpolated(xyz[:,0], xyz[:,1], xyz[:,2])
