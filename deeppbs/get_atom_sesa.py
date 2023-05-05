# standard modules
import os

# third party modules
import numpy as np

# geobind modules
from .run_msms import runMSMS
from .structure_data import StructureData

def getAtomSESA(atoms, prefix='area', mesh=None, clean=True, hydrogens=True, quiet=True, feature_name='sesa', binary=False, threshold=0.1, **kwargs):
    if isinstance(atoms, StructureData):
        atoms = atoms.atom_list
    
    if mesh is None:
        # run MSMS
        af = runMSMS(atoms, prefix, '.', area_only=True, hydrogens=hydrogens, quiet=quiet, clean=clean, **kwargs)
        
        # read in area file
        SE = open(af)
        SE.readline()
        count = 0
        for i in range(len(atoms)):
            if (not hydrogens) and atoms[i].element == 'H':
                continue
            sesa = float(SE.readline().strip().split()[1])
            if binary:
                sesa = int(sesa > threshold)
            atoms[i].xtra[feature_name] = sesa
            count += 1
        
        # clean up
        if clean:
            os.remove(af)
    else:
        from .map_point_features_to_structure import mapPointFeaturesToStructure
        
        # map vertex areas to atoms
        vertex_areas = np.zeros(len(mesh.vertices))
        np.add.at(vertex_areas, mesh.faces[:, 0], mesh.area_faces/3)
        np.add.at(vertex_areas, mesh.faces[:, 1], mesh.area_faces/3)
        np.add.at(vertex_areas, mesh.faces[:, 2], mesh.area_faces/3)
        
        mapPointFeaturesToStructure(mesh.vertices, atoms,  vertex_areas, feature_name)
        
    return feature_name

