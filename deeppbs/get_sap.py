
# geobind packages
from .data import data
from .get_atom_kdtree import getAtomKDTree

def getSAP(structure, distance=5.0, ns=None,
         bonds=None, residue_hydrophobicity=None, standard_area=None, side_chain_atoms=None,
         area_key='sesa', impute_hydrogens=False, feature_name='sap'
    ):
    # Reads in the file pdbid-protein.pdb and computes the SAP score for
    # all standard residues on the protein surface. Non-standard should 
    # be assigned a score of none. Non-standard residue atoms are ignored
    # and are not considered in SAP calculations.
    # Arguments:
    # pdbid:       structure name/identifier
    # distance:    distance cut-off for neighbor search
    #-------------------------------------------------------------------
    
    # check for data we need
    if bonds is None:
        bonds = data.covalent_bond_data
    
    if residue_hydrophobicity is None:
        residue_hydrophobicity = data.residue_hydrophobicity
    
    if standard_area is None:
        if area_key == 'sesa':
            standard_area = data.standard_sesa
        elif area_key == 'sasa':
            standard_area = data.standard_sasa
        else:
            raise ValueError("Must supply `standard_area` if `area_key` is not 'sesa' or 'sasa'!")
    
    if side_chain_atoms is None:
        side_chain_atoms = data.chem_components
    
    # check if we're given a KDTree object
    if ns is None:
        ns = structure.atom_KDTree
    
    # compute the atom SAP scores
    for chain in structure.get_chains():
        for residue in chain:
            for a in residue:
                if impute_hydrogens and a.element == 'H':
                    continue
                if a.xtra[area_key] <= 0.0:
                    a.xtra[feature_name] = 0.0
                    continue
                center = a.get_coord()
                sap = 0.0
                neighbors = ns.search(center, distance, level='A')
                for n in neighbors:
                    if impute_hydrogens and n.element == 'H':
                        continue
                    if n.xtra[area_key] <= 0.0:
                        continue
                    nname = n.get_name().strip()
                    nres = n.get_parent()
                    nresn = nres.get_resname().strip()
                    # select side-chain atoms only
                    if nname in side_chain_atoms[nresn]['side_chain_atoms']:
                        if standard_area[nresn][nname] == 0.0:
                            sap += residue_hydrophobicity[nresn] # assume a RASA of 1.0
                        else:
                            sap += residue_hydrophobicity[nresn]*min(1.5, n.xtra[area_key]/standard_area[nresn][nname])
                a.xtra[feature_name] = sap
    
    if impute_hydrogens:
        # use parent atom as hydrogen SAP value if we exluded them
        for chain in structure.get_chains():
            for residue in chain:
                resn = residue.get_resname().strip()
                for a in residue:
                    if a.element != 'H':
                        continue
                    aname = a.get_name().strip()
                    parent_atom = bonds[resn][aname]['bonded_atoms'][0]
                    if parent_atom in residue:
                        a.xtra[feature_name] = residue[parent_atom].xtra[feature_name]
    
    return feature_name

