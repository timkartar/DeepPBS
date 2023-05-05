# third party packages
import numpy as np

# geobind packages
from .get_residue_id import getResidueID
from .get_atom_kdtree import getAtomKDTree
from .data import data

def getCV(structure, radius, residue_ids=None, ns=None, feature_name="cv", impute_hydrogens=False, bonds=None):
    """get CV values for every residue in residues"""
    
    if residue_ids is None:
        # use all residues in the structure
        residue_ids = []
        for chain in structure.get_chains():
            for residue in chain:
                residue_ids.append(getResidueID(residue))
    
    if ns is None:
        # create a KDtree for structure atoms
        ns = structure.atom_KDTree
    
    for resID in residue_ids:
        cid, num, ins = resID.split('.')
        residue = structure.get_residue((' ', int(num), ins), cid)
        for atom in residue:
            if impute_hydrogens and atom.element == 'H':
                continue
            vector = np.zeros(3)
            count = 0
            neighbors = ns.search(atom.get_coord(), radius, level='A')
            for n in neighbors:
                if n == atom:
                    continue
                if impute_hydrogens and n.element == 'H':
                    continue
                dx = atom.get_coord() - n.get_coord()
                vector += dx/(np.linalg.norm(dx) + 1e-5)
                count += 1
            atom.xtra[feature_name] = 1 - np.linalg.norm(vector)/(count + 1e-5)
    
    if impute_hydrogens:
        # use parent atom as hydrogen CV value if we exluded them
        if bonds is None:
            # use default bond data
            bonds = data.covalent_bond_data
        
        for resID in residue_ids:
            cid, num, ins = resID.split('.')
            residue = structure.get_residue((' ', int(num), ins), cid)
            resn = residue.get_resname().strip()
            for atom in residue:
                if atom.element != 'H':
                    continue
                aname = atom.get_name().strip()
                parent_atom = bonds[resn][aname]['bonded_atoms'][0]
                if parent_atom in residue:
                    atom.xtra[feature_name] = residue[parent_atom].xtra[feature_name]
    
    return feature_name
