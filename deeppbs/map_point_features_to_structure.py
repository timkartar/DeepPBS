import numpy as np

from .one_hot_encode import oneHotEncode
from .get_atom_kdtree import getAtomKDTree

def mapPointFeaturesToStructure(points, atom_list, features, feature_names, kdtree=None, reduce_method='sum', impute=False, impute_value=0.0):
    if features.ndim == 1:
        features = features.reshape(-1, 1)
    
    if isinstance(feature_names, str):
        feature_names = [feature_names]
    
    # build a KDTree for structure if one is not provided
    if kdtree is None:
        kdtree = getAtomKDTree(atom_list)
    
    # reset values if atom list was used prior
    for atom in atom_list:
        for fn in feature_names:
            atom.xtra[fn] = []
    
    # find nearest-neighbor atoms and accumulate features
    dist, ind = kdtree.query(points)
    F = len(feature_names)
    for i in range(len(ind)):
        ai = ind[i]
        pi = i
        for j in range(F):
            fn = feature_names[j]
            fv = features[pi, j]
            atom_list[ai].xtra[fn].append(fv)
    
    # reduce
    if reduce_method == 'sum':
        reduce_fn = np.sum
    elif reduce_method == 'mean':
        reduce_fn = np.mean
    elif reduce_method == 'max':
        reduce_fn = np.max
    
    for atom in atom_list:
        for fn in feature_names:
            if len(atom.xtra[fn]) == 0:
                if impute:
                    atom.xtra[fn] = impute_value
                else:
                    del atom.xtra[fn]
            else:
                atom.xtra[fn] = reduce_fn(atom.xtra[fn])

def mapVertexProbabilitiesToStructure(vertices, atom_list, P, nc, level='A', kdtree=None, vertex_weights=None, reduce_method='mean'):
    
    # one hot encode if given labels vector
    if P.ndim == 1:
        P = oneHotEncode(P, nc)
    
    # build a KDTree for structure if one is not provided
    if kdtree is None:
        kdtree = getAtomKDTree(atom_list)
    
    # reset values if atom list was used prior
    for atom in atom_list:
        atom.xtra['p'] = np.zeros(nc)
        atom.xtra['v_weight'] = 0
    
    # find nearest-neighbor atoms and add probabilities
    if vertex_weights is None:
        vertex_weights = np.ones(len(vertices))
    
    dist, ind = kdtree.query(vertices)
    for i in range(len(ind)):
        ai = ind[i]
        vi = i
        atom_list[ai].xtra['p'] += P[vi]*vertex_weights[vi]
        atom_list[ai].xtra['v_weight'] += vertex_weights[vi]
    
    # normalize probabilities
    for atom in atom_list:
        if atom.xtra['v_weight'] > 0:
            atom.xtra['p'] = atom.xtra['p']/atom.xtra['v_weight']
    
    # return list of entities with aggregated probabilities
    if level == 'A':
        atom_dict = {atom.get_full_id(): atom.xtra['p'] for atom in atom_list}
        return atom_dict
    elif level == 'R':
        residue_dict = {}
        # aggregate over atom probabilities
        for atom in atom_list:
            residue = atom.get_parent()
            residue_id = residue.get_full_id()
            if residue_id not in residue_dict:
                residue_dict[residue_id] = []
            residue_dict[residue_id].append(atom.xtra['p'])
        
        # reduce residue probabilities
        for residue_id in residue_dict:
            residue_dict[residue_id] = np.stack(residue_dict[residue_id])
            if reduce_method == 'mean':
                residue_dict[residue_id] = np.mean(residue_dict[residue_id], axis=0)
            elif reduce_method == 'max':
                residue_dict[residue_id] = np.max(residue_dict[residue_id], axis=0)
        
        return residue_dict

