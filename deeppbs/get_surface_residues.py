# geobind modules
from .get_residue_id import getResidueID

def getSurfaceResidues(structure, area_key='sesa', threshold=0.0, hydrogens=False):
    atoms = structure.get_atoms()
    
    # get residue area values
    residueArea = {} # dict to map rid to an area value
    for atom in atoms:
        if (not hydrogens) and atom.element == 'H':
            continue

        rid = getResidueID(atom.get_parent())
        if rid not in residueArea:
            residueArea[rid] = 0.0
        residueArea[rid] += atom.xtra[area_key]
    
    # determine surface residues
    surface_residues = []
    for rid in residueArea:
        if residueArea[rid] > threshold:
            surface_residues.append(rid)
    
    return surface_residues
