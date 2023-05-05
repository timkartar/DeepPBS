# builtin modules
import subprocess
import os

# geobind modules
from .strip_hydrogens import stripHydrogens
from .structure_data import StructureData
from .data import data

def getAtomChargeRadius(structure, prefix='structure', hydrogens=True, keepPQR=False, source="AMBER", min_radius=0.6):
    """ Assign atomic parameters to every atom in a protein chain. Values are stored in atom.xtra 
    dictionary. The following keys are used
    ------------------------------------------------------------------------------------------------
    radius - van der Waals radius
    charge - atomic effective charge
    
    ARGUMENTS:
        structure - a geobind Structure object
    """
    if source == "AMBER":
        for atom in structure.get_atoms():
            residue = atom.get_parent()
            resn = residue.get_resname().strip()
            atmn = atom.name.strip()
            
            if resn in data.standard_residues:
                # protein checks
                if resn == "HIS":
                    # check protonation state
                    if "HD1" in residue and "HE2" in residue:
                        resn = "HIP"
                    elif "HD1" in residue:
                        resn = "HID"
                    elif "HE2" in residue:
                        resn = "HIE"
                    else:
                        resn = "HIP"
                
                # check for C-terminus
                if "OXT" in residue:
                    resn = "C"+resn
                
                # check for N-terminus
                if ("H1" in residue) or ("H2" in residue) or ("H3" in residue):
                    resn = "N"+resn
            
            if resn in data.standard_DNA_nucleotides:
                # DNA checks
                
                # check for 5' end
                if "HO5'" in residue:
                    resn = resn+"5"
                
                # check for 3' end
                if "HO3'" in residue:
                    resn = resn+"3"
                
                # additional atom naming conventions
                if atmn == "H5'":
                    atmn = "H5'1"
                
                if atmn == "H5''":
                    atmn = "H5'2"
                    
                if atmn == "H2'":
                    atmn = "H2'1"
                    
                if atmn == "H2''":
                    atmn = "H2'2"
                    
                if atmn == "HO5'":
                    atmn = "H5T"
                
                if atmn == "HO3'":
                    atmn = "H3T"
                
                if atmn == "OP1":
                    atmn = "O1P"
                    
                if atmn == "OP2":
                    atmn = "O2P"
                
                if atmn == "OP3":
                    atmn = "O2P" # substitute for existing atom
                    
            atom.xtra["radius"] = data.AMBER[resn][atmn]["radius"]
            atom.xtra["charge"] = data.AMBER[resn][atmn]["charge"]
    elif source == "freesasa":
        # assign radius only
        from. get_atom_sasa import Radius
        R = Radius()
        
        for atom in structure.get_atoms():
            residue = atom.get_parent()
            resn = residue.get_resname().strip()
            atmn = atom.name.strip()
            
            atom.xtra["radius"] = R.radius(resn, atmn)
    elif source == "EDTSurf":
        # assign radius only
        
        for atom in structure.get_atoms():
            atmn = atom.name.strip()
            elem = atom.element
            
            if atmn in data.vdw_radii_edtsurf:
                atom.xtra["radius"] = data.vdw_radii_edtsurf[atmn]
            elif elem in data.vdw_radii_edtsurf:
                atom.xtra["radius"] = data.vdw_radii_edtsurf[elem]
            else:
                atom.xtra["radius"] = data.vdw_radii_edtsurf['*']

