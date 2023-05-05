from Bio.PDB import NeighborSearch
import numpy as np
import sys

def countContacts(protein, pdb, V_dna, dna_mask):
    try:
        cid = pdb.split("_")[1].split(".")[0]
        ns = NeighborSearch(list(protein[cid].get_atoms()))
    except:
        ns = NeighborSearch(list(protein.get_atoms()))
    V_dna = V_dna[dna_mask,:,:].reshape(-1,3)
    count = 0
    for item in V_dna:
        out = ns.search(item, 5)
        count += len(out)

    return np.array([count])
