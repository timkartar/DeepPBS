import os

def getDSSP(model, dssp_map=None, feature_name='secondary_structure', formatstr="{}({})", clean=True):
    try:
        from Bio.PDB.DSSP import DSSP
    except ModuleNotFoundError:
        raise ModuleNotFoundError("The module 'Bio.PDB.DSSP' is required for this functionality!")
    
    if dssp_map is None:
        # map eight ss types to three
        dssp_map = {
            "H": formatstr.format(feature_name, "H"),
            "G": formatstr.format(feature_name, "H"),
            "I": formatstr.format(feature_name, "H"),
            "E": formatstr.format(feature_name, "S"),
            "B": formatstr.format(feature_name, "L"),
            "T": formatstr.format(feature_name, "L"),
            "S": formatstr.format(feature_name, "L"),
            "-": formatstr.format(feature_name, "L")
        }
    
    # Write a PDB file
    pdb_file = model.save()
    
    # run DSSP using the DSSP class from BioPython
    dssp = DSSP(model, pdb_file)
    
    # store secondary structure in each atom property dict
    keys = list(sorted(set(dssp_map.values())))
    for chain in model:
        cid = chain.get_id()
        for residue in chain:
            rid = residue.get_id()
            dkey = (cid, rid)
            
            if dkey in dssp:
                ss = dssp_map[dssp[dkey][2]]
            else:
                ss = dssp_map['-']
            
            for atom in residue:
                atom.xtra[keys[0]] = 0.0
                atom.xtra[keys[1]] = 0.0
                atom.xtra[keys[2]] = 0.0
                atom.xtra[ss] = 1.0
    
    if clean:
        os.remove(pdb_file)
    
    return keys

