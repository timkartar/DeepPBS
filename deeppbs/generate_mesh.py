# builtin modules
import os
import logging

# third party modules
import numpy as np

# geobind modules
from .run_nanoshaper import runNanoShaper
from .run_msms import runMSMS
from .structure import StructureData 

def generateMesh(structure, 
        prefix=None, basedir=None, clean=True, hydrogens=True, quiet=True, 
        method='nanoshaper', selection=None, fallback=None, **kwargs
    ):
    def _makeMesh(method, structure):
        if method == 'nanoshaper':
            # Run NanoShaper
            mesh = runNanoShaper(structure.atom_list, prefix, basedir, clean=clean, hydrogens=hydrogens, quiet=quiet, **kwargs)
        elif method == 'msms':
            # Run MSMS
            mesh = runMSMS(structure.atom_list, prefix, basedir, clean=clean, hydrogens=hydrogens, quiet=quiet, **kwargs)
        elif method == 'edtsurf':
            # Run EDTSurf
            pdbfile = structure.save('tmp.pdb')
            mesh = runEDTSurf(pdbfile, prefix, basedir, clean=clean, quiet=quiet, **kwargs)
        
        return mesh
    
    # Check what we have been given
    if isinstance(structure, str):
        # Check if PQR file exists
        if not os.path.exists(structure):
            raise ValueError("Can not find PQR file: {}".format(structure))
            
        # Load a PQR file
        if prefix is None:
            prefix = ".".join(os.path.basename(structure).split('.')[:-1]) # strip the file extension
        structure = StructureData(structure, name=prefix)
        if entity_id:
            # Choose which part of the structure we want to use to generate a mesh
            structure = structure.slice(structure, selection)
        
        # add charge/radius
        for atom in structure.atom_list:
            atom.xtra["radius"] = atom.bfactor
            atom.xtra["charge"] = atom.occupancy
    else:
        # assume this is a StructureData object
        level = structure.get_level()
        if prefix is None:
            # Use the structure id as the prefix
            if level == 'S':
                prefix = structure.get_id()
            elif level == 'M':
                prefix = structure.get_parent().get_id()
            elif level == 'C':
                prefix = structure.get_parent().get_parent().get_id()
    
    # decide where mesh is saved
    if basedir == '.' or basedir is None:
        basedir = os.getcwd()
    
    try:
        mesh = _makeMesh(method, structure)
    except FileNotFoundError as e:
        if method == fallback:
            raise e
        logging.info("%s failed to produce any output. Falling back to %s", method, fallback)
        mesh = _makeMesh(fallback, structure)
    
    return mesh
