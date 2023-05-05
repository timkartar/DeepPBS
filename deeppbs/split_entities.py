# geobind modules
from .data import data as D

def splitEntities(structure, regexes=None, atom_mapper=None, mi=0):
    """Docstring"""
    if regexes is None:
        regexes = D.regexes
    
    pro = []
    lig = []
    for chain in structure[mi]:
        for residue in chain:
            resname = residue.get_resname().strip()
            if regexes['PROTEIN']['STANDARD_RESIDUES'].search(resname):
                pro.append(residue.get_full_id())
            else:
                if atom_mapper is not None:
                    if atom_mapper.testResidue(residue):
                        lig.append(residue.get_full_id())
                else:
                    lig.append(residue.get_full_id())
    
    return structure.slice(structure, pro, name='{}_protein'.format(structure.name)), structure.slice(structure, lig, name='{}_ligand'.format(structure.name))
