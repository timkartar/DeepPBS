# built in modules
import os
import shutil
import subprocess
from random import choice
from string import ascii_letters

# third party modules
import numpy as np

# geobind modules
from .strip_hydrogens import stripHydrogens
from .structure_data import StructureData

def padCoordinates(pqrFile):
    padded = open("tmp.pqr", "w")
    with open(pqrFile) as FH:
        for line in FH:
            if line[0:4] == "ATOM":
                s = line[:30]
                e = line[54:]
                x = line[30:38].strip()
                y = line[38:46].strip()
                z = line[46:54].strip()
                padded.write("{}{:>9s}{:>9s}{:>9s}{}".format(s, x, y, z, e))
            else:
                padded.write(line)
    padded.close()
    shutil.move("tmp.pqr", pqrFile)

def tempFileName(prefix, ext):
    return "%s.%s" % (prefix + ''.join(choice(ascii_letters) for i in range(15)), ext)

def runPDB2PQR(structure, replace_hydrogens=False, command="pdb2pqr", add_charge_radius=True, keep_pqr=True, min_radius=0.6, structure_name=None, pad_coords=True):
    try:
        from Bio.PDB import PDBParser
    except ModuleNotFoundError:
        raise ModuleNotFoundError("The dependency 'BioPython' is required!")
    
    if structure_name is None:
        structure_name = structure.name
    
    FNULL = open(os.devnull, 'w')
    # check if PDB2PQR is installed
    rc = subprocess.call(['which', command], stdout=FNULL, stderr=FNULL)
    if rc:
        raise Exception("Command {} not found when trying to run PDB2PQR!".format(command))
    
    if replace_hydrogens:
        # strip any existing hydrogens - add new ones with PDB2PQR
        stripHydrogens(structure)
    
    # Write chain to temp file
    tmpFile = tempFileName(structure_name, 'pdb')
    pqrFile = "{}.pqr".format(structure_name)
    structure.save(tmpFile)
    
    # Run PDB2PQR
    if add_charge_radius:
        chain_flag = "--keep-chain"
    else:
        chain_flag = ""
    subprocess.call([
            command,
            '--ff=AMBER',
            chain_flag,
            tmpFile,
            pqrFile
        ],
        stdout=FNULL,
        stderr=FNULL
    )
    FNULL.close()
    
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    if not os.path.exists(pqrFile):
        raise FileNotFoundError("No PQR file was produced ({}). Try manually running PDB2PQR on the pbdfile file '{}' and verify output.".format(pqrFile, tmpFile))
    structure = parser.get_structure(structure_name, pqrFile)
    model = structure[0]
    
    # Get radius and charge from PQR file
    for line in open(pqrFile):
        if line[0:4] != "ATOM":
            continue
        cid = line[21]
        num = int(line[22:26].strip())
        ins = line[26]
        rid = (" ", num, ins)
        crg = float(line[55:62].strip())
        vdw = float(line[63:69].strip())
        atm = line[12:16].strip()
        if vdw == 0.0:
            vdw = min_radius # 0 radius atoms causes issues - set to a minimum of 0.6
        if rid in model[cid]:
            if add_charge_radius and (atm in model[cid][rid]):
                model[cid][rid][atm].xtra["charge"] = crg 
                model[cid][rid][atm].xtra["radius"] = vdw
    structure = StructureData(model, name=structure_name)
    
    # clean up
    os.remove(tmpFile)
    if keep_pqr:
        if pad_coords:
            padCoordinates(pqrFile)
        return structure, pqrFile
    else:
        os.remove(pqrFile)
        return structure

