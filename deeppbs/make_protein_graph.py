from Bio.PDB import PDBParser, NeighborSearch
import numpy as np
import sys
import json
import os

atom_keys = ['C','O','N','S']

BOND_LENGTH_THRESHOLD = 1.96 # threshold determined based on https://www.umass.edu/microbio/rasmol/rasbonds.htm

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_data/bonds.json")

def makeProteinGraph(model, feature_names=None, skip_hydrogens=True):
    #global atom_map
    V = []
    X = []
    E = []
    
    if skip_hydrogens:
        atoms = list(filter(lambda a: a.element != 'H', model.get_atoms()))
    else:
        atoms = list(model.get_atoms())
    
    ns = NeighborSearch(atoms)

    bonds = json.load(open(DATA_PATH,"r"))   

    vectors = []
    for atom in atoms:
        atom_feature = np.zeros(len(atom_keys))
        idx = atom_keys.index(atom.element)
        atom_feature[idx] = 1
        
        # get covalently bonded neighbors
        res = atom.get_parent()
        res_bonds = bonds[res.get_resname()]

        bonded_atoms = []
        for i in res_bonds[atom.get_name()]:
            if i == "OXT":
                bonded_atoms = ns.search(center=atom.coord, radius=BOND_LENGTH_THRESHOLD, level="A")[1:] 
            elif i[0] == "H":
                continue
            else:
                try:
                    bonded_atoms.append(res[i])
                except:
                    print("Missing protein atom", i)
                
        vs = np.array([(atom.coord - bonded_atom.coord) for bonded_atom in bonded_atoms])
        v = np.mean(vs, axis=0)
        v = v/np.linalg.norm(v)
        vectors.append(v)
        
        for bonded_atom in bonded_atoms:
            E.append([atoms.index(atom), atoms.index(bonded_atom)]) # this has O(N^2) complexity... should be improved. Can replace with hash table for O(N)

        V.append(atom.coord)
        
        x = list(atom_feature) + [atom.xtra[f] for f in feature_names] # list of features
        X.append(x)
    
    vectors = np.array(vectors)
    
    return np.array(V), np.array(X), np.array(E), vectors
