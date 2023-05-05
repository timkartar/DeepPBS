def stripHydrogens(structure):
    """Strip all hydrogen atoms from the given model.
    
    Parameters
    ----------
    model: BioPython MODEL object
        The model to be stripped of hydrogens.
    """
    if structure.get_level() == 'R':
        rm = []
        for atom in structure:
            if atom.element == 'H':
                rm.append(atom.get_id())
        for aid in rm:
            structure.detach_child(aid)
    else:
        for residue in structure.get_residues():
            rm = []
            for atom in residue:
                if atom.element == 'H':
                    rm.append(atom.get_id())
            for aid in rm:
                residue.detach_child(aid)
