import numpy as np

def getAtomKDTree(atoms, engine='scipy'):
    # check if we're passed a list of atom objects or a higher level entity
    if not isinstance(atoms, list):
        atoms = [atom for atom in atoms.get_atoms()]
    
    # check which KDTree we want to construct. These have non-equilvalent methods
    if engine.lower() == 'scipy':
        from scipy.spatial import cKDTree
        coords = np.array([atom.coord for atom in atoms])
        
        return cKDTree(coords)
    elif engine.lower() == 'biopython':
        from Bio.PDB import NeighborSearch
        
        return NeighborSearch(atoms)
    else:
        raise ValueError("Unknown engine for generating atom KDTree: {}".format(engine))
