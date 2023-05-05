
def getAtomDepth(model, vertices=None, feature_name="depth"):
    try:
        from Bio.PDB.ResidueDepth import min_dist, get_surface
    except ModuleNotFoundError:
        raise ModuleNotFoundError("The module 'Bio.PDB.ResidueDepth' is required for this functionality!")
    
    if vertices is None:
        vertices = get_surface(model.structure)
    
    for atom in model.get_atoms():
        coord = atom.coord
        dist = min_dist(vertices, coord)
        
        atom.xtra[feature_name] = dist
    
    return feature_name
