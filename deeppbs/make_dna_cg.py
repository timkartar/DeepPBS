import numpy as np
from .dna_encodings import seqToOneHot
from .data import data 
import sys

# atom positions in order from near backbone to near helix axis
MAJOR_GROOVE_ATOMS = {
    "DA": ["N7", "N6"],
    "DC": ["N4", "C5"],
    "DG": ["N7", "O6"],
    "DT": ["O4", "C5"]
}

MINOR_GROOVE_ATOMS = {
    "DA": ["N3", "C2"],
    "DC": ["O2"],
    "DG": ["N3", "C2"],
    "DT": ["O2"]
}

#TODO: better to hard-code these explicity
SUGAR_ATOMS = []

fn_label_dict = { 
        ('DA', 'N7'): 0, # acceptor
        ('DA', 'N6'): 1, # donor
        ('DT', 'O4'): 0,
        ('DT', 'C5'): 2, # methyl
        ('DA', 'N3'): 0,
        ('DA', 'C2'): 3, # non-polar
        ('DT', 'O2'): 0,
        ('DG', 'N7'): 0,
        ('DG', 'O6'): 0,
        ('DC', 'C5'): 3,
        ('DC', 'N4'): 1,
        ('DG', 'N3'): 0,
        ('DG', 'C2'): 1,
        ('DC', 'O2'): 0
        }
major_dict = {
        'DA': [0,1,0,2],
        'DC': [3,1,0,0],
        'DG': [0,0,1,3],
        'DT': [2,0,1,0]
}

minor_dict = {
	'DA': [0,3,0],
        'DC': [0,1,0],
        'DG': [0,1,0],
        'DT': [0,3,0]
}

def get_fn(base1):

    fn = np.array([[0.0]*6]*11)
    fn[[0,10],4] = 1
    fn[[1,9],5] = 1
    

    for i in range(4):
        fn[i+2, major_dict[base1][i]] = 1
    
    for i in range(3):
        fn[i+6, minor_dict[base1][i]] = 1
    return fn

def fix_resname(rname):
    if rname == 'DI':
        return 'DG'
    elif rname == 'DU':
        return 'DT'
    else:
        return rname
def getP(res):
    if "P" in res and "OP1" in res and "OP2" in res:
        ### return P vector ###
        v1 = res["OP1"].coord - res["P"].coord
        v2 = res["OP2"].coord - res["P"].coord
        v = (v1 + v2)
        v = v/np.linalg.norm(v)

        return res["P"].coord, v
    
    return None, [0,0,0]

def getS(res):
    coords = []
    for atom in res.get_atoms():
        if "'" in atom.get_name():
            coords.append(atom.coord)
    
    try:
        v = res["C5'"].coord - res["C4'"].coord
    except:
        v1 = res["C4'"].coord - res["C3'"].coord
        v2 = res["C4'"].coord - res["O4'"].coord
        v = np.cross(v1, v2)
    v = v/np.linalg.norm(v)

    return np.mean(coords, axis=0).tolist(), v

def get_groove(res1, res2, groove="major"): ## unused, is not symmetric
    if groove == "major":
        groove_dict = MAJOR_GROOVE_ATOMS
    elif groove == "minor":
        groove_dict = MINOR_GROOVE_ATOMS
    
    ret = []
    fn = []
    for anchor in groove_dict[res1.get_resname()]:
        fn_toappend = [0.0]*6
        for atom in res1:
            if atom.get_name() == anchor:
                ret.append(atom.coord.tolist())
                fn_toappend[fn_label_dict[(res1.get_resname(),anchor)]] = 1
                fn.append(fn_toappend)
                break
    
    for anchor in groove_dict[res2.get_resname()][::-1]:
        fn_toappend = [0.0]*6
        for atom in res2:
            if atom.get_name() == anchor:
                ret.append(atom.coord.tolist())
                fn_toappend[fn_label_dict[(res2.get_resname(),anchor)]] = 1
                fn.append(fn_toappend)
                break
    return ret, fn

def get_groove_symmetric(res1, res2):
    C11 = None
    C12 = None

    N1 = None
    N9 = None
    

    res1_atom_names = [atom.get_name() for atom in res1.get_atoms()]
    res2_atom_names = [atom.get_name() for atom in res2.get_atoms()]
    
    if 'N9' in res1_atom_names:
        purine = res1
        pyrimidine = res2
    else:
        purine = res2
        pyrimidine = res1
        

    for atom in res1.get_atoms():
        if "C1'" in atom.get_name():
            C11 = atom.coord
        if 'N9' in res1_atom_names:
            if "N9" in atom.get_name():
                N9 = atom.coord
        else:
            if "N1" in atom.get_name():
                N1 = atom.coord

   
    
    for atom in res2.get_atoms():
        if "C1'" in atom.get_name():
            C12 = atom.coord
        if 'N9' in res2_atom_names:
            if "N9" in atom.get_name():
                N9 = atom.coord
        else:
            if "N1" in atom.get_name():
                N1 = atom.coord

    try: #DEBUG
        point_on_plane = (N1 + N9)/2
    except:
        if N1 == None:
            point_on_plane = N9
        elif N9 == None:
            point_on_plane = N1
    
    v = point_on_plane - C11
    u = C12 - C11

    v_proj_u = np.dot(v, u)/(np.linalg.norm(u)**2)

    perp = v - v_proj_u*u

    perp_dir = perp/np.linalg.norm(perp)

    major_groove_start = C11 + 3.745*perp_dir - (1.54)*u/np.linalg.norm(u) ## estimated shift by major groove positions mean
    major_groove_u = u + 2*(1.54)*u/np.linalg.norm(u)
    major_groove = []
    for i in range(4):
        major_groove.append(major_groove_start + major_groove_u*(i+1)/5)

    ## add 1A vertical shift for middle two beads
    major_groove[1] += perp_dir
    major_groove[2] += perp_dir

    minor_groove = []
    for i in range(3):
        minor_groove.append(C11 + u*(i+1)/4)
    
    
    ####### compute vectors ###########

    cg = np.mean(major_groove + minor_groove, axis=0)

    major_vectors = []
    minor_vectors = []

    for item in major_groove:
        vector = item - cg
        vector = vector/np.linalg.norm(vector)
        major_vectors.append(vector)
    
    for item in minor_groove:
        vector = item - cg
        vector = vector/np.linalg.norm(vector)
        minor_vectors.append(vector)


    return major_groove, minor_groove, major_vectors, minor_vectors

def fill_5primephosphate(cg_dna,  vectors, missing_phosphates):
    for item in missing_phosphates:
        if item == 0:
            S1 = cg_dna[item,1]
            S2 = cg_dna[item+1,1]
            P = cg_dna[item+1,0]
            vectors[item,0] = vectors[item+1,0]
            cg_dna[item,0] = P + S1 - S2

        if item == (len(cg_dna) - 1):
            S1 = cg_dna[item,9]
            S2 = cg_dna[item - 1, 9]
            P = cg_dna[item - 1,10]
            vectors[item,10] = vectors[item-1,0]
            
            cg_dna[item,10] = P + S1 - S2

    return cg_dna, vectors


def makeDNACG(dna, helix, mi=0, symmetric=True):
    cg_dna = []
    seq = []
    fn = []
    vectors = []
    #p_fn = [0.0,0.0,0.0,0.0,1,0.0]
    #s_fn = [0.0,0.0,0.0,0.0,0.0,1]
    
    missing_phosphates = []
    missing_sugars = []
    for i in range(helix["length"]):
        id1 = helix["ids1"][i].split(".")
        id2 = helix["ids2"][i].split(".")
        
        base1 = dna[mi][id1[0]][int(id1[1])]
        base2 = dna[mi][id2[0]][int(id2[1])]
        
        seq.append([fix_resname(base1.get_resname())[-1], fix_resname(base2.get_resname())[-1]])
        
        P1, Pv1 = getP(base1)
        P2, Pv2 = getP(base2)
        
        if P1 is None:
            P1 = [0.0, 0.0, 0.0]
            missing_phosphates.append(i)
        if P2 is None:
            P2 = [0.0, 0.0, 0.0]
            missing_phosphates.append(i)
        
        S1, Sv1 = getS(base1)
        S2, Sv2 = getS(base2)
        
        if S1 is None:
            S1 = [0.0, 0.0, 0.0]
            missing_sugars.append(i)
        if S2 is None:
            S2 = [0.0, 0.0, 0.0] 
            missing_sugars.append(i)
        
        
        
        if symmetric:
            major_groove_positions, minor_groove_positions, major_vectors, minor_vectors = get_groove_symmetric(base1, base2)
        
        vectors.append([Pv1, Sv1] + major_vectors + minor_vectors + [Sv2, Pv2])
        cg_dna.append([P1, S1] + major_groove_positions + minor_groove_positions + [S2, P2])
        fn.append(get_fn(fix_resname(base1.get_resname())))
    

    cg_dna = np.array(cg_dna)
    vectors = np.array(vectors)

    cg_dna, vectors = fill_5primephosphate(cg_dna, vectors, missing_phosphates)   
        
    fn = np.array(fn)

        
    seq = np.array(seq).T.tolist()
    seq[0] = seqToOneHot("".join(seq[0]))
    seq[1] = seqToOneHot("".join(seq[1][::-1])) ## make second strand seq 5' -> 3'

    return cg_dna, np.array(seq), fn, vectors

