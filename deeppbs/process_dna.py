# built-in modules
import json
import re
import os
import subprocess
import math
import collections
from itertools import permutations

# third-party utils
import networkx as nx
from networkx.algorithms import bipartite
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import NeighborSearch
from Bio.PDB.PDBIO import Select
from Bio.PDB import PDBIO

import numpy as np
from sklearn import linear_model

# geobind_dna modules
from .data import data as D
import warnings
warnings.filterwarnings('ignore')

# Used to determine if a DNA helical axis is linear, curved in plane or
# curved out-of-plane
__RANK_DISCRIMINATOR = 0.12

# Various cut-off values
__STACKING_OVERLAP_CUTOFF = 0.5
__INTERBASE_ANGLE_CUTOFF = 45.0
__LINK_DISTANCE_CUTOFF = 1.85

# Name maping
nucNameMap = D.long_to_short
bpNameMap = {
    "WC": "watson-crick",
    "Hoogsteen": "hoogsteen",
    "wobble": "wobble"
}

# A map for pair rules. Since different pair types have different pair
# rules, this allows for classificaiton of mismatches based on pair
# type.
pairRules = {
    "watson-crick": {
        "A": ["T", "t", "U", "u"],
        "C": ["G", "g"],
        "G": ["C", "c"],
        "T": ["A", "a"],
        "U": ["A", "a"],
        "a": ["T", "t"],
        "c": ["G", "g"],
        "g": ["C", "c"],
        "t": ["A", "a"],
        "u": ["A", "a"]
    },
    "hoogsteen": {
        "A": ["T", "t", "U", "u"],
        "C": ["G", "g"],
        "G": ["C", "c"],
        "T": ["A", "a"],
        "U": ["A", "a"],
        "a": ["T", "t"],
        "c": ["G", "g"],
        "g": ["C", "c"],
        "t": ["A", "a"],
        "u": ["A", "a"]
    }
}

namedChemicalModifications = {
    "5CM": '5_methylated_cytosine'
}

# Regexes
HBOND_ATOM_RE = re.compile('([A-Z]+[0-9]*[\'\*]*)[(\[]?')
HBOND_DIST_RE = re.compile('\[([0-9]\.[0-9]+)\]')

class EntitySelect(Select):
    def __init__(self, ids):
        self.ids = ids
    
    def accept_residue(self, residue):
        return getID(residue=residue) in self.ids

def getID(*args, sep_char='.', **kwargs):
    if 'residue' in kwargs:
        ch = kwargs['residue'].get_parent().get_id()
        rid = kwargs['residue'].get_id()
        fields = [ch, str(rid[1]), rid[2]]
        return '.'.join(fields)
    else:
        args = [str(_) for _ in args]
        return sep_char.join(args)

def getHash(*args, sep_char=None):
    args = list(args)
    # Check that no tuples are present
    for item in args:
        if isinstance(item, tuple):
            for t in item:
                args.append(t)
            args.remove(item)
    if sep_char:
        return sep_char.join(sorted(args))
    else:
        return tuple(sorted(args))

def runDSSR(structure, quiet=True):
    """Run DSSR on given PDB file.
    
    Parameters
    ----------
    prefix: string
        The file prefix.
    """
    
    prefix = 'dna'
    if not isinstance(structure, str):
        file_name = "{}.tmp.pdb".format(prefix)
        structure.save(file_name)
    else:
        file_name = structure
    
    args = ["x3dna-dssr", "--i={}".format(file_name), "--o={}-dssr.json".format(prefix), "--json", "--more", "--idstr=long", "--non-pair"]
    if quiet:
        FNULL = open(os.devnull, 'w')
        subprocess.call(args, stdout=FNULL, stderr=FNULL)
        subprocess.call(["x3dna-dssr", "--cleanup"],stdout=FNULL, stderr=FNULL)
        FNULL.close()
    else:
        subprocess.call(args)
        subprocess.call(["x3dna-dssr", "--cleanup"])
    
    with open("{}-dssr.json".format(prefix)) as FH:
        DSSR = json.load(FH, object_pairs_hook=collections.OrderedDict)
    
    return DSSR

def getNucleotideById(model, nid, hetflag=" "):
    """Retrieve a Biopython nucleotide object from given 'model' by the 
    nucleotide id 'nid'.
    
    Parameters
    ----------
    model: Model object 
        A Biopython Model object which stores a DNA structure.
    nid: string
        The nucleotide identifier.
    """
    ch, num, ins = nid.split('.')
    rid = (hetflag, int(num), ins)
    
    return model[ch][rid]

def getNucleotideId(nucleotide):
    """Returns a nucleotide id string given a Biopython nucleotide 
    object.
    
    Parameters
    ----------
    nucleotide: Nucleotide object 
        A Biopython Nucleotide object.
    """
    chain = nucleotide.get_parent().get_id()
    nid = nucleotide.get_id()
    number = str(nid[1])
    ins = nid[2]
    
    return '.'.join((chain, number, ins))

def convertId(id_string):
    """Converts a DSSR id string to standard DNAproDB format.
    
    Parameters
    ----------
    id_string: string 
        A DSSR nucleotide id string.
    """
    components = id_string.split('.')
    
    if components[5] == '':
        components[5] = ' '
    
    return '.'.join((components[2], components[4], components[5]))

def addBackboneLinkages(G, model, nuc_dict, COMPONENTS):
    """Adds edges between nodes in G, where the edge represents a 
    phosphodiester linkage.
    
    Parameters
    ----------
    G: Networkx Graph object
        Stores edges between nodes. Nucleotides are nodes and edges are
        backbone links.
    model: Biopython Model object
        The DNA structure object.
    """
    # Store all the P and Sugar atoms
    atom_list = []
    for chain in model:
        for nucleotide in chain:
            nid = getID(residue=nucleotide)
            if(nid not in nuc_dict):
                continue
            nucn = nucleotide.get_resname().strip()
            SA = COMPONENTS[nucn]["sugar_atoms_re"]
            PA = COMPONENTS[nucn]["phosphate_atoms_re"]
            for atom in nucleotide:
                aname = atom.get_name()
                if(re.search(SA, aname)):
                    atom_list.append(atom)
                elif(re.search(PA, aname)):
                    atom_list.append(atom)
    ns = NeighborSearch(atom_list)
    
    # Find nucleotides which are linked
    linkages = ns.search_all(__LINK_DISTANCE_CUTOFF, level='A')
    link_list = []
    seen = [] # store phosphate atom ids that we have already seen
    for link in linkages:
        parent0 = link[0].get_parent()
        parent1 = link[1].get_parent()
        
        # Check if this is same residue or not
        if(parent0.get_id() == parent1.get_id()):
            continue
            
        # Look for phosphorus atom
        if(link[0].element == "P"):
            Patom = link[0]
        elif(link[1].element == "P"):
            Patom = link[1]
        else:
            continue
        
        # Check if we've done this link already
        if(Patom.get_full_id() in seen):
            continue
        
        # Treat P atom as separate entitiy - find O3' and O5' atoms linked to it
        P_neighbors = ns.search(Patom.get_coord(), __LINK_DISTANCE_CUTOFF, level='A')
        P3atom = None
        P5atom = None
        for n in P_neighbors:
            if(n.element == "P" or n.element == "C"):
                continue
            name = n.get_name().strip()
            if(re.search("OP[123]|O[123]P", name)):
                continue
                
            # Look for *3' and *5' atom based on name
            if(re.search("[A-Z]+3'$", name)):
                P3atom = n
            elif(re.search("[A-Z]+5'$", name)):
                P5atom = n
            elif(P5atom is None and (n.get_parent().get_id() == Patom.get_parent().get_id())):
                # this neighbor belongs to the same residue as the phosphorus
                P5atom = n
            elif(P3atom is None and (n.get_parent().get_id() != Patom.get_parent().get_id())):
                # this neighbor belongs to a different residue
                P3atom = n
        
        if((P3atom is None) or (P5atom is None)):
            # We have failed - skip this link
            continue
        
        dist = P3atom-Patom
        Pid = Patom.get_full_id()
        P3id = getID(residue=P3atom.get_parent())
        P3atom = P3atom.get_name()
        P5id = getID(residue=Patom.get_parent())
        P5atom = Patom.get_name()
        if(P3id[0] != P5id[0]):
            # by definition, nucleotides from different chains should not be linked
            continue
        if(P3id == P5id):
            # can't have self-linkages
            continue
        
        seen.append(Pid)
        G.add_edge(P3id, P5id, type="link")
        link_list.append({
            "3p_nuc_id": P3id,
            "3p_atom": P3atom,
            "5p_nuc_id": P5id,
            "5p_atom": P5atom,
            "link_distance": float(dist)
        })
    return link_list

def getNucleotideData(nt, model, COMPONENTS):
    """Returns relevant data from a nucleotide dictionary which is 
    retrieved from the DSSR JSON output.
    
    Parameters
    ----------
    nt: dictionary
        Nucleotide dictionary from DSSR JSON output.
    """
    ins = nt['nt_id'].split('.')[5]
    if ins == '':
        ins = ' '
    
    nid = convertId(nt['nt_id'])
    data = {
        "name": nt["nt_name"].strip(),
        "name_short": nucNameMap.get(nt["nt_name"].strip(), "X"),
        "chain": nt["chain_name"],
        "number": nt["nt_resnum"],
        "ins_code": ins,
        "id": nid,
        "glycosidic_conformation": nt["baseSugar_conf"],
        "origin": nt["frame"]["origin"],
        "secondary_structure": "other",
        "modified": False,
        "chemical_name": COMPONENTS[nt["nt_name"].strip()]['_chem_comp.name']
    }
    if "is_modified" in nt:
        data["modified"] = True
    
    nucleotide = getNucleotideById(model, nid)
    
    # check if phosphate linking atom present
    ppAtom = COMPONENTS[data["name"]]["sugar-phosphate_bond"]["phosphate_atom"]
    if ppAtom in nucleotide:
        data["phosphate_present"] = True
    else:
        data["phosphate_present"] = False
    
    return data

def getPairData(pair):
    """Returns relevant data from a base pair dictionary which is 
    retrieved from the DSSR JSON output.
    
    Parameters
    ----------
    pair: dictionary
        Base pair dictionary from DSSR JSON output.
    """
    data = {
        "id1": convertId(pair["nt1"]),
        "id2": convertId(pair["nt2"]),
        "name1": pair["nt1"].split('.')[3],
        "name2": pair["nt2"].split('.')[3],
        "origin": pair["frame"]["origin"],
        "x_axis": pair["frame"]["x_axis"],
        "y_axis": pair["frame"]["y_axis"],
        "z_axis": pair["frame"]["z_axis"],
        "pair_type": bpNameMap.get(pair["name"], "other"),
        "dssr_description": pair["DSSR"],
        "hbonds": []
    }
    data["id"] = getHash(data["id1"], data["id2"])
    
    # Add hydrogen bonds
    hbs = pair["hbonds_desc"].split(",")
    for hb in hbs:
        hb1, hb2 = re.split('[-*?]', hb)
        m1 = HBOND_ATOM_RE.search(hb1)
        m2 = HBOND_ATOM_RE.search(hb2)
        atm1 = m1.group(1)
        atm2 = m2.group(1)
        m = HBOND_DIST_RE.search(hb)
        dist = m.group(1)
        data["hbonds"].append({
            "atom1": atm1,
            "atom2": atm2,
            "distance": float(dist)
        })
    
    # determine if this pair counts as a mismatch
    nt1, nt2 = re.split("[-+]", pair["bp"])
    if data["pair_type"] in pairRules:
        if nt2 in pairRules[data["pair_type"]][nt1]:
            mm = False
        else:
            mm = True
    else:
        if nt2 in pairRules["watson-crick"][nt1]:
            mm = False
        else:
            mm = True
    data["mismatched"] = mm
    
    return data

def _pointsOrientation(p1, p2, p3):
    s = (p2[1]-p1[1])*(p3[0]-p2[0]) - (p3[1]-p2[1])*(p2[0]-p1[0])
    if(s == 0):
        return 0
    elif(s > 0):
        return 1
    else:
        return 2

def _onSegment(p1, p2, p3):
    if(
        p2[0] <= max(p1[0], p3[0]) and p2[0] >= min(p1[0], p3[0]) and
        p2[1] <= max(p1[1], p3[1]) and p2[1] >= min(p1[1], p3[1])
    ):
       return True
    
    return False

def _segmentsIntersect(p1, q1, p2, q2):
    """Determine if two line segements defined by the pairs (p1, q1) and
    (p2, q2) intersect, while checking for special cases of colinearity.
    
    Parameters
    ----------
    p1: iterable, float
    q1: iterable, float
    p2: iterable, float
    q2: iterable, float
    """
    # Find the four orientations needed for general and
    # special cases
    o1 = _pointsOrientation(p1, q1, p2)
    o2 = _pointsOrientation(p1, q1, q2)
    o3 = _pointsOrientation(p2, q2, p1)
    o4 = _pointsOrientation(p2, q2, q1)
    
    # General Case
    if(o1 != o2 and o3 != o4):
        return True
    
    # Special Cases
    # p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if(o1 == 0 and _onSegment(p1, p2, q1)):
        return True
    # p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if(o2 == 0 and _onSegment(p1, q2, q1)):
        return True
    # p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if(o3 == 0 and _onSegment(p2, p1, q2)):
        return True
    # p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if(o4 == 0 and _onSegment(p2, q1, q2)):
        return True
    
    # Doesn't fall in any of the above cases
    return False

def _weight(P, a, b, c):
    """ Returns the minimum distance from the point P to the line 
    defined by a*y + b*x + c = 0
    
    Parameters
    ----------
    P: iterable, float
        Any datatype which can be accessed as P[0] and P[1] to retrieve 
        the x and y coordinates of the point.
    a: float
        y coefficient
    b: float
        x coefficient
    c: float
        y intercept
    """
    if(b is None):
        return abs(P[0] - c)
    else:
        return abs(a*P[1] + b*P[0] + c)/np.sqrt(a**2 + b**2)

def _assignNucleotideAtoms(nuc, nucid, comp, elements, edge, sg, wg, center, origin, x_axis, y_axis):
    for atom in nuc:
        atmn = atom.get_name()
        if(re.search(comp["sugar_atoms_re"], atmn) or re.search(comp["phosphate_atoms_re"], atmn)):
            # ignore sugar and phosphate atoms
            continue
        if(atmn in sg or atmn in wg):
            # skip if already included in sg or wg
            continue
        if(atmn in elements):
            if(elements[atmn] == 'H'):
                # ignore hydrogens
                continue
        else:
            if(atmn[0] == 'H'):
                # ignore hydrogens
                continue
        atmid = "{}.{}".format(nucid, atmn)
        
        # find the closest neighbor in edge
        x = np.dot(atom.coord-origin, x_axis)
        y = np.dot(atom.coord-origin, y_axis)
        point = (x, y)
        # Check number of times the point-center line segement crosses
        # the minor groove edge
        crossing_count = 0
        for i in range(len(edge)-1):
            if(_segmentsIntersect(point, center, edge[i], edge[i+1])):
                crossing_count += 1
        
        if(crossing_count % 2 == 0):
            sg.append(atmn)
        else:
            wg.append(atmn)

def _addAtomEdges(G, nuc, nid, comp, elements, origin, x_axis, y_axis, a, b, c):
    """Adds edges between atoms based on bonds described in
    '_chem_comp_bond' and computes weights for each edge based on the 
    distance of the center point of each edge from the line described 
    by ax +by + c = 0
    
    Parameters
    ----------
    G: networkx Graph object
        describes the atom edges of the base-pair
    nuc: Biopython Residue object
        nucleotide object
    nid: string
        nucleotide ID string
    comp: dict
        PDB component entry corresponding to the nucleotide
    elements: dict
        lookup dictionary describing the element of each atom in 'nuc'
    origin: numpy array 
        origin of the base-pair reference frame
    x_axis: numpy array
        x-axis of the base-pair reference frame
    y_axis: numpy array
        y-axis of the base-pair reference frame
    a: float
        x-coefficient of a line in the minor groove region
    b: float
        y-coefficient of a line in the minor groove region
    c: float
        constant of a line in the minor groove region
    """
    for i in range(len(comp["_chem_comp_bond.atom_id_1"])):
        atm1 = comp["_chem_comp_bond.atom_id_1"][i]
        atm2 = comp["_chem_comp_bond.atom_id_2"][i]
        if(re.search(comp["sugar_atoms_re"], atm1) or re.search(comp["phosphate_atoms_re"], atm1)):
            # skip sugar and phosphate atom edges
            continue
        elif(re.search(comp["sugar_atoms_re"], atm2) or re.search(comp["phosphate_atoms_re"], atm2)):
            # skip sugar and phosphate atom edges
            continue
        elif(elements[atm1] == 'H' or elements[atm2] == 'H'):
            # skip hydrogens
            continue
        elif(atm1 in nuc and atm2 in nuc):
            x1 = np.dot(nuc[atm1].coord-origin, x_axis)
            y1 = np.dot(nuc[atm1].coord-origin, y_axis)
            x2 = np.dot(nuc[atm2].coord-origin, x_axis)
            y2 = np.dot(nuc[atm2].coord-origin, y_axis)
            weight = _weight(((x1+x2)/2, (y1+y2)/2), a, b, c)
            G.add_edge(
                "{}.{}".format(nid,atm1),
                "{}.{}".format(nid,atm2),
                weight=weight**4
            )

def partitionGrooveAtoms(model, pair, nuc_dict, COMPONENTS):
    """Given a helical base pair, partitions the base atoms into minor 
    and major groove atoms for each nucleotide.
    
    Parameters
    ----------
    model: Biopython Model object
        model object containing the base pair
    pair: dict
        pair JSON dictionary
    nuc_dict: dict
        nucleotide lookup, keyed by 'nuc_id'
    COMPONENTS: dict
        dictionary of PDB component descriptions
     """
    n1 = getNucleotideById(model, pair["id1"])
    n2 = getNucleotideById(model, pair["id2"])
    nid1 = pair["id1"]
    nid2 = pair["id2"]
    
    comp1 = COMPONENTS[pair["name1"]]
    comp2 = COMPONENTS[pair["name2"]]
    y_axis = np.array(pair["y_axis"])
    x_axis = np.array(pair["x_axis"])
    origin = np.array(pair["origin"])
    
    start = comp1["sugar-base_bond"]["base_atom"]
    end = comp2["sugar-base_bond"]["base_atom"]
    
    # compute center point as a point perpendicular to the line joining
    # the projected N1(9) and N9(1) atom coordinates of the two paired 
    # bases, opposite of the direction of the x-axis unit vector.
    NS = np.array([
        np.dot(n1[start].coord-origin, x_axis),
        np.dot(n1[start].coord-origin, y_axis)
    ])
    NE = np.array([
        np.dot(n2[end].coord-origin, x_axis),
        np.dot(n2[end].coord-origin, y_axis)
    ])
    NN = NE-NS
    
    # Direction perpendicular to line NN
    center = np.array([-NN[1]/np.sqrt(np.dot(NN,NN)), NN[0]/np.sqrt(np.dot(NN,NN))])
    
    # Determine direction to move center point based on various edge 
    # cases
    coef = 8 # distance to place groove center
    if(pair["dssr_description"][2] == "+" and pair["dssr_description"][1] == pair["dssr_description"][3]):
        # Reflect if pair is in syn-anti conformation
        coef *= -1
    if(center[0] > 0):
        # Reflect if in the MG direction
        coef *= -1
    center *= coef
    a = 1
    if(NN[0] > 0.001):
        b = -NN[1]/NN[0]
        c = -(center[1] + b*center[0])
    else:
        b = None
        c = center[0]
    
    A = nx.Graph()
    # Add element assignments to each atom in components
    elements1 = {}
    for i in range(len(comp1["_chem_comp_atom.type_symbol"])):
        elements1[comp1["_chem_comp_atom.atom_id"][i]] =  comp1["_chem_comp_atom.type_symbol"][i]
    
    elements2 = {}
    for i in range(len(comp2["_chem_comp_atom.type_symbol"])):
        elements2[comp2["_chem_comp_atom.atom_id"][i]] =  comp2["_chem_comp_atom.type_symbol"][i]
    
    _addAtomEdges(A, n1, nid1, comp1, elements1, origin, x_axis, y_axis, a, b, c)
    _addAtomEdges(A, n2, nid2, comp2, elements2, origin, x_axis, y_axis, a, b, c)
    
    for hb in pair["hbonds"]:
        atm1 = hb["atom1"]
        atm2 = hb["atom2"]
        if re.search(comp1["sugar_atoms_re"], atm1) or re.search(comp1["phosphate_atoms_re"], atm1):
            continue
        elif re.search(comp2["sugar_atoms_re"], atm2) or re.search(comp2["phosphate_atoms_re"], atm2):
            continue
        x1 = np.dot(n1[atm1].coord-origin, x_axis)
        y1 = np.dot(n1[atm1].coord-origin, y_axis)
        x2 = np.dot(n2[atm2].coord-origin, x_axis)
        y2 = np.dot(n2[atm2].coord-origin, y_axis)
        weight = _weight(((x1+x2)/2, (y1+y2)/2), a, b, c)
        A.add_edge(
            "{}.{}".format(nid1, atm1),
            "{}.{}".format(nid2, atm2),
            weight=weight**4
        )
    
    # Get the shortest path
    start = "{}.{}".format(nid1, start)
    end = "{}.{}".format(nid2, end)
    try:
        p = nx.shortest_path(A, source=start, target=end, weight='weight')
    except:
        # nucleotide bases are disjoint - add a fake hydrogen bond
        mindist = 99999
        for a1 in n1:
            aname1 = a1.name
            if(re.search(comp1["sugar_atoms_re"], aname1) or re.search(comp1["phosphate_atoms_re"], aname1)):
                # skip sugar and phosphate atom edges
                continue
            for a2 in n2:
                aname2 = a2.name
                if(re.search(comp2["sugar_atoms_re"], aname2) or re.search(comp2["phosphate_atoms_re"], aname2)):
                    # skip sugar and phosphate atom edges
                    continue
                dist = a2-a1
                if(dist < mindist):
                    mindist = dist
                    atmpair = (a1, a2)
        atm1 = atmpair[0]
        atm2 = atmpair[1]
        x1 = np.dot(atm1.coord-origin, x_axis)
        y1 = np.dot(atm1.coord-origin, y_axis)
        x2 = np.dot(atm2.coord-origin, x_axis)
        y2 = np.dot(atm2.coord-origin, y_axis)
        weight = _weight(((x1+x2)/2, (y1+y2)/2), a, b, c)
        A.add_edge(
            "{}.{}".format(nid1, atm1.name),
            "{}.{}".format(nid2, atm2.name),
            weight=weight**4
        )
        p = nx.shortest_path(A, source=start, target=end, weight='weight')
    
    # add projected coordinate of each node in mG edge
    edge_points = []
    edge_points.append(NS - 100*NN)
    for node in p:
        f = node.split('.')
        nid = '.'.join(f[0:3])
        atm = f[3]
        if(nid == nid1):
            coord = n1[atm].coord
        else:
            coord = n2[atm].coord
        x = np.dot(coord-origin, x_axis)
        y = np.dot(coord-origin, y_axis)
        edge_points.append((x, y))
    edge_points.append(NE + 100*NN)
    
    # store assignments for each atom
    sg = {
        nid1: [],
        nid2: []
    }
    wg = {
        nid1: [p.pop(-0).split('.')[-1]],
        nid2: [p.pop(-1).split('.')[-1]]
    }
    
    # add atoms from edge to minor groove
    for atom in p:
        f = atom.split('.')
        nid = '.'.join(f[0:3])
        atm = f[3]
        sg[nid].append(atm)
    
    # assign rest of base atoms
    _assignNucleotideAtoms(n1, nid1, comp1, elements1, edge_points, sg[nid1], wg[nid1], center, origin, x_axis, y_axis)
    _assignNucleotideAtoms(n2, nid2, comp2, elements2, edge_points, sg[nid2], wg[nid2], center, origin, x_axis, y_axis)
    
    # remove N1/N9 atoms - don't classify these
    wg[nid1].remove(wg[nid1][0])
    wg[nid2].remove(wg[nid2][0])
    
    nuc_dict[nid1]["groove_atoms"] = {
        'sg': sg[nid1],
        'wg': wg[nid1]
    }
    nuc_dict[nid2]["groove_atoms"] = {
        'sg': sg[nid2],
        'wg': wg[nid2]
    }

def getStackData(stack):
    """Returns relevant data from a stacking interaction dictionary 
    which is retrieved from the DSSR JSON output.
    
    Parameters
    ----------
    stack: dictionary
        Stacking interaction dictionary from DSSR JSON output.
    """
    data = {
        "id1": convertId(stack["nt1"]),
        "id2": convertId(stack["nt2"]),
        "name1": stack["nt1"].split('.')[2],
        "name2": stack["nt2"].split('.')[2],
        "mindist": stack["min_baseDist"]
    }
    data["id"] = getHash(data["id1"], data["id2"])
    
    return data

def addPairs(G, dssr):
    """Iterates through the DSSR JSON output and adds each base-pair as 
    an edge in G.
    
    Parameters
    ----------
    G: Networkx Graph object
        Stores edges between nodes. Nucleotides are nodes and edges are
        base-pairings.
    dssr: dictionary
        The parsed DSSR JSON output. 
    """
    pair_list = []
    if "pairs" in dssr:
        for pair in dssr["pairs"]:
            #if(pair["interBase_angle"] > __INTERBASE_ANGLE_CUTOFF):
            #    continue
            G.add_edge(convertId(pair["nt1"]), convertId(pair["nt2"]), type="pair")
            pair_list.append(getPairData(pair))
    
    return pair_list

def addStacks(G, dssr):
    """Iterates through the DSSR JSON output and adds each stack as 
    an edge in G.
    
    Parameters
    ----------
    G: Networkx Graph object
        Stores edges between nodes. Nucleotides are nodes and edges are
        nucleotide stackings.
    dssr: dictionary
        The parsed DSSR JSON output. 
    """
    stack_list = []
    if("nonPairs" in dssr):
        for nonPair in dssr["nonPairs"]:
            if("stacking" in nonPair):
                if(nonPair["stacking"]["oArea"] < __STACKING_OVERLAP_CUTOFF):
                    continue
                G.add_edge(convertId(nonPair["nt1"]), convertId(nonPair["nt2"]), type="stack")
                stack_list.append(getStackData(nonPair))
    
    return stack_list

def helixSegments(sG, pG, dssr, nucExplicitNumberMap, file_name, pDict, nucleotides):
    """Tests if the stacking and pairing graphs describe a double-helix
    conformation. A double helix must have the topology of a ladder
                                N--N   -- Base Pairing
                                |  |    | Stacking
                                N--N
                                |  |
                                N--N
    but some amount of deviation from this ideal topology should be 
    allowed for. This function computes a numerical score indicating how
    close to an ideal double helix the DNA conformation is, with 1.0 
    being perfect.
    
    The score is determined by the ratio of base-pairs to cannonical 
    base pairs. A cannonical base pair is a base pair which is involved 
    in a perfect ladder-like stacking relationship.
    
    Also runs 3DNA on each helical region, as detected by DSSR, and 
    determines a helical axis for each region.
    
    Parameters
    ----------
    sG: Networkx Graph object
        Stores edges between nodes. Nucleotides are nodes and edges are
        nucleotide stacks.
    pG: Networkx Graph object
        Stores edges between nodes. Nucleotides are nodes and edges are
        base-pairings.
    dssr: dictionary
        The parsed DSSR JSON output for a DNA entity.
    nucExplicitNumberMap: dictionary
        Maps nucleotide id to sequential index.
    prefix: string
        Structure prefix string.
    pDict: dictionary
        Maps pair id to a pair JSON object.
    nucleotides: dictionary
        Maps nucleotide id to a corresponding nucleotide object.
    """
    
    segments = []
    # Parse helices found by DSSR and computes a helix score for each.
    if "helices" in dssr:
        for helix in dssr["helices"]:
            if helix["num_pairs"] < 4:
                continue
            summary = {
                "score": None,
                "multiplets": [],
                "ids1": [],
                "ids2": [],
                "sequence1": helix["strand1"],
                "sequence2": helix["strand2"],
                "length": helix["num_pairs"],
                "GC_content": None,
                "score": None,
                "helix_id": [],
                "contains_mismatches": False,
                "contains_non-wc_pairs": False,
                "contains_hoogsteen_pairs": False,
                "contains_non-standard_pairs": False,
                "chemical_modifications": {
                    "5_methylated_cytosine": False,
                    "non-standard_nucleotides": False
                }
            }
            if "helical_radius" in helix:
                summary["mean_radius"] = helix["helical_radius"]
            
            nnums = {}
            for p in helix["pairs"]:
                id1 = convertId(p["nt1"])
                id2 = convertId(p["nt2"])
                ch1, num1, _ = id1.split(".")
                ch2, num2, _ = id2.split(".")
                if ch1 not in nnums:
                    nnums[ch1] = []
                if ch2 not in nnums:
                    nnums[ch2] = []
                nnums[ch1].append(int(num1))
                nnums[ch2].append(int(num2))
                summary["ids1"].append(id1)
                summary["ids2"].append(id2)
                nucleotides[id1]["secondary_structure"] = "helical"
                nucleotides[id2]["secondary_structure"] = "helical"
                pid = getHash(id1, id2)
                if pDict[pid]["mismatched"]:
                    summary["contains_mismatches"] = True
                if pDict[pid]["pair_type"] != "watson-crick":
                    summary["contains_non-wc_pairs"] = True
                if pDict[pid]["pair_type"] == "hoogsteen":
                    summary["contains_hoogsteen_pairs"] = True
                if pDict[pid]["pair_type"] == "other":
                    summary["contains_non-standard_pairs"] = True
                
            # Set up helix id
            for ch in nnums:
                nummin = np.min(nnums[ch])
                nummax = np.max(nnums[ch])
                summary["helix_id"].append("{}{}/{}".format(ch,nummin,nummax))
            summary["helix_id"] = ".".join(summary["helix_id"])
            
            pairs = []
            for u in summary["ids1"]:
                nu = nx.node_connected_component(pG,u)
                pairs.append(list(nu))
                
            # Check pairing/stacking relationships
            cpCount = 0.0 # canonical pairs
            for p in pairs:
                if(len(p) > 2):
                    summary["multiplets"].append({
                        "ids": p,
                        "order": len(p)
                    })
                    continue
                u = p[0]
                v = p[1]
                if sG.has_node(u):
                    snu = sG.neighbors(u) # stacking neighbors of u
                else:
                    continue
                if sG.has_node(v):
                    snv = sG.neighbors(v) # stacking neighbors of v
                else:
                    continue
                apu = set([frozenset(e) for e in pG.edges(snu)]) # adjacent pairs of u
                apv = set([frozenset(e) for e in pG.edges(snv)]) # adjacent pairs of v
                
                if (len(apu) + len(apv)) > 0:
                    cpCount += 2.0*len(apu & apv)/(len(apu) + len(apv))
            
            summary["score"] = cpCount/len(pairs)
            print("Helix score: {}".format(summary["score"]))
            
            # Add shape parameters and helical axis 
            summary["shape_parameters"] = {}
            inp, out = makeInputFile(file_name, helix["pairs"], nucExplicitNumberMap)
            runAnalyze(summary, inp, out)
            mgw, MGW = runCurves(file_name, helix["pairs"], nucExplicitNumberMap)
            
            summary["shape_parameters"]["minor_groove_curves"] = mgw
            summary["shape_parameters"]["major_groove_curves"] = MGW
            
            # Compute helical axis
            summary["helical_axis"] = {}
            #getHelicalAxis(helix["pairs"], pDict, summary["helical_axis"])
            
            # Compute A-tracts
            atractRe = re.compile('(A|T(?!A)){4,}')
            Atracts = False
            if atractRe.search(summary["sequence1"]) or atractRe.search(summary["sequence2"]):
                Atracts = True
            summary["contains_A-tracts"] = Atracts
            
            # Compute GC content
            GC_content = (
                summary["sequence1"].count('G') + 
                summary["sequence1"].count('C') +
                summary["sequence2"].count('G') +
                summary["sequence2"].count('C')
            )/(2.0*summary["length"])
            summary["GC_content"] = GC_content
            
            # Update chemical modifications
            updateChemicalModifications(summary["ids1"]+summary["ids2"], nucleotides, summary["chemical_modifications"])
            
            segments.append(summary)
    
    return segments

def getHelicalAxis(pairs, pDict, aDict):
    """Computes the helical axis for a helix segment. 
    
    Parameters
    ----------
    pairs: array
        Array of JSON base-pair objects. Only used to obtain the correct
        pair IDs in the correct order.
    pDict: dictionary
        A dictionary that maps a pair ID to a pair object. These objects
        are needed because they contain the base-pair origins needed to
        compute the axis.
    aDict: dictionary
        Dictionary to store results of the fitting and other helix 
        related data.
    """
    # Store the origin coordinates
    x = []
    y = []
    z = []
    npts = len(pairs)
    for p in pairs:
        id1 = convertId(p["nt1"])
        id2 = convertId(p["nt2"])
        pid = getHash(id1, id2)
        origin = pDict[pid]["origin"]
        x.append(origin[0])
        y.append(origin[1])
        z.append(origin[2])
    
    # Determine axis geometry
    if(npts < 8):
        # Short helices always linear. 
        aDict["axis_curvature"] = 'linear'
        deg = 1
    else:
        haxis = np.array([x, y, z], dtype=float).T
        U, d, V = np.linalg.svd(haxis-haxis.mean(axis=0), full_matrices=True)
        d /= d[0]
        if(d[0] > __RANK_DISCRIMINATOR >= d[1] >= d[2]):
            aDict["axis_curvature"] = 'linear'
            deg = 1
        elif(d[0] >= d[1] > __RANK_DISCRIMINATOR >= d[2]):
            aDict["axis_curvature"] = 'in-plane'
            deg = 8
        elif(d[0] >= d[1] >= d[2] > __RANK_DISCRIMINATOR):
            aDict["axis_curvature"] = 'non-planar'
            deg = 12
        else:
            aDict["axis_curvature"] = "error"
            return
    
    # Compute a polynomial fit to the base-pair origins
    s = np.zeros(npts, dtype=float)
    for j in range(1, npts):
        dx = x[j] - x[j-1]
        dy = y[j] - y[j-1]
        dz = z[j] - z[j-1]
        s[j] = s[j-1] + math.sqrt(dx**2 + dy**2 + dz**2)
    aDict["axis_length"] = s[npts-1]
    s -= s[npts-1]/2 # center the distance measure
    si = np.linspace(s[0], s[-1], 5*npts)
    Si = []
    for n in range(1, deg+1):
        Si.append(si**n)
    Si = np.array(Si).T
    
    xi, xcoef = fitPolynomial(s, x, deg, Si)
    yi, ycoef = fitPolynomial(s, y, deg, Si)
    zi, zcoef = fitPolynomial(s, z, deg, Si)
    aDict["x_coef"] = xcoef
    aDict["y_coef"] = ycoef
    aDict["z_coef"] = zcoef
    aDict["bp_origin_x"] = x
    aDict["bp_origin_y"] = y
    aDict["bp_origin_z"] = z

def fitPolynomial(t, y, degree, ti):
    """Fit a polynomial of degree 'degree' to the data t,y and return
    a prediction for the points 'ti'.
    
    Parameters
    ----------
    t: ndarray
        Parametric coordinate for the datapoints.
    y: array
        Response variable values.
    degree: int
        Degree of the polynomial fit.
    ti: ndarray
        Numpy array of shape n,d where n is the number of datapoints
        and d is the degree of the fit.
    """
    y = np.array(y)
    T = []
    for n in range(1, degree+1):
        T.append(t**n)
    T = np.array(T).T
    
    if degree == 1:
        lambd = [0.0]
    elif degree > 1 and degree < 10:
        lambd = [0.25, 0.5, 1.0, 2.0]
    else:
        lambd = [0.25, 0.5, 1.0, 2.0, 5.0]
    
    # Try different values of lambda
    models = []
    scores = []
    for l in lambd:
        P = linear_model.Ridge(alpha=l, fit_intercept=True, normalize=False)
        P.fit(T,y)
        yp = P.predict(T)
        MSE = np.mean((y-yp)**2)
        models.append(P)
        scores.append(MSE)
        if MSE < 1.0:
            break
    
    i = np.argmin(scores)
    P = models[i]
    return P.predict(ti), [P.intercept_]+P.coef_.tolist()

def makeInputFile(file_name, pairs, nucExplicitNumberMap):
    """Generates the input file for analyze.
    
    Parameters
    ----------
    prefix: string
        The input filename prefix.
    pairs: array
        An array of JSON objects describing the base-pairing.
    nucExplicitNumberMap: dictionary
        A mapping from nucleotide id to sequential index.
    """
    prefix = file_name.split('.')[0]
    
    FH = open("{}.inp".format(prefix), "w")
    
    FH.write("{}\n".format(file_name))
    FH.write("{}.out\n".format(prefix))
    FH.write("{:>5d}\n".format(2))
    FH.write("{:>5d}\n".format(len(pairs)))
    FH.write("{:>5d} {:>5d}\n".format(1,1))
    for p in pairs:
        FH.write("{:>5d} {:>5d} {:>3d}\n".format(nucExplicitNumberMap[convertId(p["nt1"])],nucExplicitNumberMap[convertId(p["nt2"])],0))
    FH.close()
    
    return "{}.inp".format(prefix), "{}.out".format(prefix)

def ssSegments(strands, pG, nucleotides):
    """For each strand, determines the number of single-stranded 
    segments and their length. This information is returned as a list 
    of dictionaries, one dictionary per strand, and each dictionary 
    describing the found segments, if any. 
    
    A segment is a contiguous stretch of linked nucleotides where no 
    nucleotide makes any pairing to another strand, and any stacking 
    interactions are strictly within the same strand. Minimum segment 
    length is 3.
    
    Parameters
    ----------
    strands: list
        A list of dictionaries, each dictionary representing a strand 
        in the entity.
    pG: Networkx Graph object
        Stores edges between nodes. Nucleotides are nodes and edges are
        base-pairings.
    nucleotides: dictionary
        Maps nucleotide id to a corresponding nucleotide object.
    """
    # For each strand, determine single-stranded segments
    segments = []
    for s in strands:
        ids = set(s["ids"])
        ssnuc = [True]*len(s["ids"]) # store whether nucleotide is single-stranded
        
        # Get pairing info
        index = 0
        start = 0
        count = 0 # number of consecutive 'S' pairs
        spseg = [] # store segments of continuous 'S' pairs
        for node in s["ids"]:
            # Check that node is not in a helical region
            if nucleotides[node]['secondary_structure'] == "helical":
                ssnuc[index] = False
                index +=1 
                continue
            
            # Check that nodes don't make pairs with other strands and 
            # check for continuous stretches of intra-strand pairs.
            if(pG.has_node(node)):
                pn = set(pG.neighbors(node))
                if(len(pn - ids) == 0 and len(pn) == 1):
                    pt = 'S'
                else:
                    pt = 'O'
                    ssnuc[index] = False
            else:
                pt = 'N'
            
            if(pt == 'S'):
                if(count == 0):
                    start = index
                count += 1
            else:
                if(count > 2): # self-pairing segment threshold (3)
                    spseg.append((start,start+count))
                count = 0
            index += 1
        if count > 2:
            spseg.append((start,start+count))
        
        # Remove self-pairing segments
        for seg in spseg:
            for i in range(seg[0], seg[1]):
                ssnuc[i] = False
        
        # Identify single stranded segments
        count = 0
        start = 0
        index = 0
        for ss in ssnuc:
            if ss:
                if count == 0:
                    start = index
                count += 1
            else:
                if count > 2:
                    segment = {
                        "ids": s["ids"][start:start+count],
                        "sequence": s["sequence"][start:start+count],
                        "length": count,
                        "chain_id": s["chain_id"]
                    }
                    segments.append(segment)
                count = 0
            index += 1
        if count > 2:
            segment = {
                "ids": s["ids"][start:start+count],
                "sequence": s["sequence"][start:start+count],
                "length": count,
                "chain_id": s["chain_id"]
            }
            segments.append(segment)
    
    # Update nucleotide secondary structure
    for seg in segments:
        for nid in seg["ids"]:
            nucleotides[nid]['secondary_structure'] = 'single-stranded'
    
    return segments

def runCurves(file_name, pairs, nucExplicitNumberMap):
    """Runs curves 5.3"""
    prefix = file_name.split('.')[0]
    
    crv_file = "{}.crv".format(prefix)
    lis_file = "{}.lis".format(prefix)
    grp_file = "{}_grp.pdb".format(prefix)
    pdb_file = "{}.pdb".format(prefix)
    tmp_file = "{}.temp".format(pdb_file)
    
    # remove Curves junk if present
    if os.access(lis_file, os.R_OK):
        os.remove(lis_file)
    if os.access(grp_file, os.R_OK):
        os.remove(grp_file)
    
    # rename DNA residues to make curves happy
    os.rename(pdb_file, tmp_file)
    TMP = open(pdb_file,'w')
    rc = subprocess.call(["sed", r"s/D\([ACGT]\) /\1  /g", tmp_file], stdout=TMP)
    TMP.close()
    
    # run Curves
    CRV = open(crv_file, "w")
    CRV.write("&inp file={0:}, comb=.t., fit=.t., grv=.t., \n     lis={1:}, pdb={1:}_grp, &end\n".format(pdb_file, prefix))
    CRV.write("2 {0:} -{0:} 0 0\n".format(len(pairs)))
    pairs1 = []
    pairs2 = []
    for p in pairs:
        pairs1.append(str(nucExplicitNumberMap[convertId(p["nt1"])]))
        pairs2.append(str(nucExplicitNumberMap[convertId(p["nt2"])]))
    CRV.write(" {}\n".format(" ".join(pairs1)))
    CRV.write(" {}\n".format(" ".join(pairs2)))
    CRV.write("0.0 0.0 0.0 0.0")
    CRV.close()
    
    CRV = open(crv_file)
    rc += subprocess.call(["Cur5"], stdin=CRV)
    CRV.close()
    os.rename(tmp_file, pdb_file)
    
    mgw = []
    MGW = []
    # Use Curves groove width values if possible
    if rc == 0 and os.access(lis_file, os.R_OK):
        LISFILE = open(lis_file).readlines()
        i = 0
        # Get to the groove width data
        sub_levels = None
        while i < len(LISFILE):
            m = re.search('Atom defining backbone: P\s+\d+\s+levels,\s+(\d+)\s+sub-levels', LISFILE[i])
            if m:
                sub_levels = int(m.group(1))
                i += 5
                break
            else:
                i += 1
        if sub_levels is None:
            print("Could not read groove widths from Curves output!")
            return mgw, MGW
        else:
            sub_levels = sub_levels // 2
        
        # Read in the groove width data
        mg = []
        MG = []
        li = []
        lc = 0
        while i < len(LISFILE):
            # Get mgw value
            if re.search('\d*\.?\d+', LISFILE[i][14:19]):
                mg.append(float(LISFILE[i][14:19].strip(' *')))
            else:
                mg.append('nan')
            # Get MGW value
            if re.search('\d*\.?\d+', LISFILE[i][42:47]):
                MG.append(float(LISFILE[i][42:47].strip(' *')))
            else:
                MG.append('nan')
            # Check if we are at a level start/end
            if re.search('[ACGT]\s+\d+\s+\d', LISFILE[i]):
                li.append(lc)
            lc += 1
            i += 1
        mg = np.array(mg, dtype=np.float32)
        MG = np.array(MG, dtype=np.float32)
        for i in li:
            s = max(0, i-sub_levels)
            e = min(len(mg), i+sub_levels+1)
            print(s, e)
            m = np.nanmean(mg[s:e])
            M = np.nanmean(MG[s:e])
            mgw.append('NA' if np.isnan(m) else float(m))
            MGW.append('NA' if np.isnan(M) else float(M))
    else:
        print("Curves 5.3 failed to produce output!")
    
    # clean up
    if os.access(lis_file, os.R_OK):
        os.remove(lis_file)
    if os.access(grp_file, os.R_OK):
        os.remove(grp_file)
    os.remove(crv_file)
    
    return mgw, MGW

def runAnalyze(helix, fileName, output, quiet=False):
    """Runs analyze on the given input file and extracts shape 
    parameters from the 3DNA output.
    
    Parameters
    ----------
    helix: dictionary
        A dictionary which stores data for a helix.
    fileName: string
        The inputfile for analyze.
    outpt: string
        Name of the ouput file analyze will produce.
    """
    NRE = '[+-]?\d*\.?\d+|-+'
    if quiet:
        FNULL = open(os.devnull, 'w')
        rc = subprocess.call(["analyze", fileName], stdout=FNULL, stderr=FNULL)
        FNULL.close()
    else:
        rc = subprocess.call(["analyze", fileName])
    
    # Extract Data from analyze output
    if rc == 0 and os.access(output, os.R_OK):
        FH = open(output)
        
        # Regular Expressions for parsing various components of output
        breakRe1 = re.compile('^\*+$')
        breakRe2 = re.compile('^\s+~+$')
        pairRe = re.compile('\s+bp\s+Shear\s+Stretch\s+Stagger\s+Buckle\s+Propeller\s+Opening')
        pairExtractRe = re.compile('\s+\d+\s+[ACGUTacgtu][+-][ACGUTacgtu]\s+({0})\s+({0})\s+({0})\s+({0})\s+({0})\s+({0})'.format(NRE))
        stepRe = re.compile('\s+step\s+Shift\s+Slide\s+Rise\s+Tilt\s+Roll\s+Twist')
        stepExtractRe = re.compile('\s+\d+\s+[ACGUTacgtu]+/[ACGUTacgtu]+\s+({0})\s+({0})\s+({0})\s+({0})\s+({0})\s+({0})'.format(NRE))
        classRe = re.compile('Structure classification:')
        formRe = re.compile('step\s+Xp\s+Yp\s+Zp\s+XpH\s+YpH\s+ZpH\s+Form')
        grooveRe = re.compile('P-P\s+Refined\s+P-P\s+Refined')
        grooveExtractRe = re.compile('\s+\d+\s+[ACGUTacgut]+/[ACGUTacgut]+\s+(?:{0})\s+({0})\s+(?:{0})\s+({0})'.format(NRE))
        classInfo = ""
        bp_param_names = ['shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening']
        bs_param_names = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
        
        # Add shape parameter arrays
        for p in bp_param_names:
            helix["shape_parameters"][p] = []
        for p in bs_param_names:
            helix["shape_parameters"][p] = []
            
        # Extract Base-Pair Parameters
        for line in FH:
            if pairRe.search(line):
                for line in FH:
                    if breakRe1.search(line) or breakRe2.search(line):
                        break
                    else:
                        mtch = pairExtractRe.search(line)
                        if mtch:
                            for i in range(len(bp_param_names)):
                                try:
                                    helix['shape_parameters'][bp_param_names[i]].append(float(mtch.group(i+1)))
                                except ValueError:
                                    helix['shape_parameters'][bp_param_names[i]].append('NA')
                        else:
                            print("Error parsing 3DNA output!", fileName)
                break
        
        # Extract Base-Pair Step Parameters
        for line in FH:
            if stepRe.search(line):
                for line in FH:
                    if breakRe1.search(line) or breakRe2.search(line):
                        break
                    else:
                        mtch = stepExtractRe.search(line)
                        if mtch:
                            for i in range(len(bs_param_names)):
                                try:
                                    helix['shape_parameters'][bs_param_names[i]].append(float(mtch.group(i+1)))
                                except ValueError:
                                    helix['shape_parameters'][bs_param_names[i]].append('NA')
                        else:
                            print("Error parsing 3DNA output!", fileName)
                break
        
        # Get Helix Classification
        structureClass = None
        helixClass = None
        for line in FH:
            if classRe.search(line):
                for line in FH:
                    if breakRe1.search(line):
                        break
                    else:
                        classInfo += line
                break
        if re.search('This is a right-handed nucleic acid structure', classInfo):
            helixClass = 'Right-Handed'
        elif re.search('This is a left-handed Z-form structure', classInfo):
            helixClass = 'Left-Handed'
            structureClass = 'Z-DNA'
        
        # Get DNA Form Info
        bp_step_classification = []
        for line in FH:
            if formRe.search(line):
                for line in FH:
                    if breakRe1.search(line):
                        break
                    else:
                        bp_step_classification.append(line[59:].strip())
                break
        if structureClass is None:
            Bfraction = bp_step_classification.count('B')/float(len(bp_step_classification))
            Afraction = bp_step_classification.count('A')/float(len(bp_step_classification))
            if Bfraction >= 0.5:
                structureClass = 'B-DNA'
            elif Afraction >= 0.5:
                structureClass = 'A-DNA'
            else:
                structureClass = 'other'
        helix['classification'] = structureClass
        helix['handedness'] = helixClass
        
        # Get Groove Widths
        helix['shape_parameters']['minor_groove_3dna'] = []
        helix['shape_parameters']['major_groove_3dna'] = []
        for line in FH:
            if grooveRe.search(line):
                for line in FH:
                    if breakRe1.search(line):
                        break
                    else:
                        mtch = grooveExtractRe.search(line)
                        try:
                            helix['shape_parameters']['minor_groove_3dna'].append(round(float(mtch.group(1))-5.8,3))
                        except ValueError:
                            helix['shape_parameters']['minor_groove_3dna'].append('NA')
                        try:
                            helix['shape_parameters']['major_groove_3dna'].append(round(float(mtch.group(2))-5.8,3))
                        except ValueError:
                            helix['shape_parameters']['major_groove_3dna'].append('NA')
                break
        os.remove(output)
        os.remove(fileName)
    else:
        print("analyze failed to run.", output)
        os.remove(fileName)

def getPairMap(data, pair_dict, nucleotides):
    pair_map = {}
    # Add non-WC pairs to pair table
    for j in range(len(nucleotides)):
        nid = nucleotides[j]
        pair_map[nid] = None
        if nid in pair_dict:
            # Add WC or Hoogsteen basepairs
            for pair in pair_dict[nid]:
                if pair["pair_type"] == "watson-crick" or pair["pair_type"] == "hoogsteen":
                    if nid == pair["id1"]:
                        pair_map[nid] = pair["id2"]
                    else:
                        pair_map[nid] = pair["id1"]
                    break
    return pair_map

def flatten(l):
    return reduce(lambda a,b: a + (flatten(b) if hasattr(b, '__iter__') else [b]), l, [])

def makePairTuples(ordered_ids, pair_map):
    pair_tuples = []
    for i in range(len(ordered_ids)):
        if(pair_map[ordered_ids[i]] == None):
            pair_tuples.append((i, None))
        else:
            pair_tuples.append((i, ordered_ids.index(pair_map[ordered_ids[i]])))
    return pair_tuples

def updateChemicalModifications(nuc_ids, nuc_dict, chem_dict):
    for nid in nuc_ids:
        if(nuc_dict[nid]["modified"]):
            chem_dict[namedChemicalModifications.get(nuc_dict[nid]["name"], "non-standard_nucleotides")] = True

def processDNA(structure, quiet=True):
    """Main function which processes a DNA structures and extracts
    various information from it. Each model in the structure is stored
    as a JSON object in an array.
    
    Parameters
    ----------
    prefix: string
        Filename prefix string of a PDB file which contains only 
        nucleotides and has all ligands, protein, and solvent removed.
    N: int
        Number of models to process.
    quiet: boolean (True)
        Set to true to suppress output from external programs.
    """
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    
    # Run DSSR and load the output json
    data = runDSSR(structure, quiet)
    model = structure[0]
    
    # Array to store DNA data. Each entry is a dict which contains a 
    # list of entities, where each entity describes a connected component
    # of the nucleotide graph.
    OUT_DATA = []
    
    # Add nucleotides
    nuc_list = []
    nuc_dict = {}
    G = nx.MultiGraph()
    for nt in data["nts"]:
        n = getNucleotideData(nt, model, D.chem_components)
        nuc_list.append(n)
        G.add_node(n["id"])
        nuc_dict[n["id"]] = n
    
    # Add edges in graphs
    links = addBackboneLinkages(G, model, nuc_dict, D.chem_components)
    pairs = addPairs(G, data)
    stacks = addStacks(G, data)
    
    # Map pair id to pair
    pair_dict = {}
    for p in pairs:
        pair_dict[p["id"]] = p
        if p["id1"] not in pair_dict:
            pair_dict[p["id1"]] = []
        if p["id2"] not in pair_dict:
            pair_dict[p["id2"]] = []
        pair_dict[p["id1"]].append(p)
        pair_dict[p["id2"]].append(p)
    
    # Process entities. Each entity is a connected subgraph of G.
    entities = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    
    # Get chain info
    all_chains = {}
    for nt in data["nts"]:
        cid = nt["chain_name"]
        nid = convertId(nt["nt_id"])
        name = nucNameMap.get(nt["nt_name"], "X")
        if not (cid in all_chains):
            all_chains[cid] = {
                "ids": [],
                "index": [],
                "sequence": []
            }
        all_chains[cid]["ids"].append(nid)
        all_chains[cid]["sequence"].append(name)
        all_chains[cid]["index"].append(int(nt["index_chain"]))
    for key in list(data["dbn"].keys())[1:]:
        cid = key[-1]
        all_chains[cid]["ids"] = [n for i,n in sorted(zip(all_chains[cid]["index"], all_chains[cid]["ids"]))]
        all_chains[cid]["sequence"] = [n for i,n in sorted(zip(all_chains[cid]["index"], all_chains[cid]["sequence"]))]
        all_chains[cid]["index"].sort()
    
    # Output data structure
    out = {
        "pairs": pairs,
        "stacks": stacks,
        "nucleotides": nuc_list,
        "links": links,
        "entities": [],
        "num_entities": None,
        "num_chains": len(data["chains"]),
        "chains": [],
        "num_nucleotides": len(nuc_list), 
        "ignored_entities": []
    }
    for chain in data["chains"]:
        cid = chain[-1]
        out["chains"].append({
            "id": cid,
            "num_nucleotides": data["chains"][chain]["num_nts"],
            "nucleotide_ids": all_chains[cid]["ids"],
            "sequence": "".join(all_chains[cid]["sequence"]),
        })
    
    # Iterate over entities
    ecount = 0
    strandNum = {}
    for entity in entities:
        if entity.number_of_nodes() < 2:
            # single nucleotide fragments are ignored.
            out["ignored_entities"].append(entity.nodes())
            continue
        eout = {
            "type": None, # helix, ssDNA, hybrid, or other
            "nucleotides": list(entity.nodes()),
            "strands": [],
            "pairs": [],
            "stacks": [],
            "helical_segments": None,
            "single-stranded_segments": None,
            "id": [],
            "num_helical_segments": None,
            "num_strands": None,
            "num_single-stranded_segments": None,
            "chemical_modifications": {
                "5_methylated_cytosine": False,
                "non-standard_nucleotides": False
            },
            "num_nucleotides": None
        }
        lG = nx.Graph() # link sub-graph
        sG = nx.Graph() # stack sub-graph
        pG = nx.Graph() # pair sub-graph
        
        # Add edges to sub-graphs
        for u,v,d in entity.edges(data=True):
            if d["type"] == "link":
                lG.add_edge(u, v)
            elif d["type"] == "pair":
                pG.add_edge(u, v)
            elif d["type"] == "stack":
                sG.add_edge(u, v)
            else:
                raise ValueError("Unknown edge type: {}".format(d["type"]))
        
        # Add pairs to entity
        for p in pG.edges():
            eout["pairs"].append(getHash(p))
        
        # Add stacks to entity
        for s in sG.edges():
            eout["stacks"].append(getHash(s))
        
        # Add strands to entity
        strands = [lG.subgraph(c).copy() for c in nx.connected_components(lG)]
        strand_dict = {}
        for s in strands:
            P3end = None
            P5end = None
            sd = {
                "ids": [],
                "sequence": "",
                "length": None,
                "chain_id": None,
                "strand_id": None,
                "5p_end": None,
                "3p_end": None,
                "GC_content": None,
                "chemical_modifications": {
                    "5_methylated_cytosine": False,
                    "non-standard_nucleotides": False
                }
            }
            # Identify the 5' and 3' ends of the strand
            for u,d in s.degree():
                if d == 1:
                    # check if 5' or 3' end
                    for l in links:
                        if l['3p_nuc_id'] == u:
                            P3end = u
                            break
                        elif l['5p_nuc_id'] == u:
                            P5end = u
                            break
            sd["5p_end"] = P5end
            sd["3p_end"] = P3end
            
            # Iterate over strand in 3'->5' order
            node = P3end
            j = 0
            seen_chains = set()
            while j <= len(s):
                sd["ids"].append(node)
                ch, num, ins = node.split('.')
                nid = (' ', int(num), ins)
                sd["sequence"] += nucNameMap.get(model[ch][nid].get_resname().strip(), 'X')
                seen_chains.add(ch)
                strand_dict[node] = sd
                if node == P5end:
                    break
                neighbors = s.neighbors(node)
                for n in neighbors:
                    if(n in sd["ids"]):
                        continue
                    else:
                        node = n
                j += 1
            sd["length"] = len(sd["ids"])
            
            # add strand and chain id
            if len(seen_chains) != 1:
                print("Invalid number of chains per strand! Strands must have exactly one chain ID.", prefix)
            cid = seen_chains.pop()
            snum = strandNum.get(cid, 1)
            sd["strand_id"] = cid+str(snum)
            sd["chain_id"] = cid
            eout["id"].append(sd["strand_id"])
            strandNum[cid] = snum + 1
            
            # Compute GC content
            sd["GC_content"] = (sd["sequence"].count('G') + sd["sequence"].count('C'))/float(sd["length"])
            
            # Update chemical modifications
            updateChemicalModifications(sd["ids"], nuc_dict, sd["chemical_modifications"])
            
            # add strand to entity
            eout["strands"].append(sd)
            strand_dict[sd["strand_id"]] = sd
        
        # Add any missing nucleotides as one-nt strands
        leftover = set(eout["nucleotides"]) - set(lG.nodes())
        for nid in leftover:
            ch, num, ins = nid.split('.')
            snum = strandNum.get(ch, 1)
            strandNum[ch] = snum + 1
            rid = (' ', int(num), ins)
            sd = {
                "ids": [nid],
                "sequence": nucNameMap.get(model[ch][rid].get_resname().strip(), 'X'),
                "length": 1,
                "chain_id": ch,
                "strand_id": ch+str(snum),
                "5p_end": nid,
                "3p_end": nid
            }
            eout["id"].append(sd["strand_id"])
            eout["strands"].append(sd)
            strand_dict[nid] = sd
            strand_dict[sd["strand_id"]] = sd
        eout["num_strands"] = len(eout["strands"])
        
        # Create PDB file for entity using Biopython selection
        io = PDBIO()
        io.set_structure(model)
        ename = "{}_entity_{}.pdb".format('dna', ecount)
        io.save(ename, EntitySelect(eout["nucleotides"]))
        
        # Map residue id to sequential residue number
        ent = parser.get_structure(ename, ename)
        nucSeqMap = {}
        sequential = 1
        for chain in ent[0]:
            for residue in chain:
                nucSeqMap[getID(residue=residue)] = sequential
                sequential += 1
        edata = runDSSR(ename)
        
        # Determine entity type
        hseg = helixSegments(sG, pG, edata, nucSeqMap, ename, pair_dict, nuc_dict)
        sseg = ssSegments(eout["strands"], pG, nuc_dict)
        ntotal = len(eout["nucleotides"])
        if len(hseg) > 1:
            # Structures with multiple helices automatically
            # classified as 'other'
            etype = "other"
        elif len(hseg) == 1:
            # Determine fraction of helical and single-stranded
            # nucleotides
            htotal = 0
            stotal = 0
            if hseg[0]["score"] <  0.6:
                htype = "irregular_helix"
            elif hseg[0]["score"] >= 0.6 and hseg[0]["score"] <= 0.9:
                htype = "imperfect_helix"
            else:
                htype = "perfect_helix"
            htotal = 2*hseg[0]["length"]
            for s in sseg:
                stotal += s["length"]
            htotal = htotal/float(ntotal)
            stotal = stotal/float(ntotal)
            if htotal > 0.6:
                etype = htype
            elif stotal > 0.6:
                etype = "ssDNA"
            elif htotal + stotal >= 0.6:
                etype = htype+"/ssDNA"
            else:
                etype = "other"
        else:
            # Determine fraction of single-stranded nucleotides
            stotal = 0
            for s in sseg:
                stotal += s["length"]
            stotal = stotal/float(ntotal)
            if stotal > 0.6:
                etype = "ssDNA"
            else:
                etype = "other"
        
        # Add groove assignments to helical segment pairs
        for helix in hseg:
            for k in range(len(helix["ids1"])):
                pid = getHash(helix["ids1"][k], helix["ids2"][k])
                partitionGrooveAtoms(model, pair_dict[pid], nuc_dict, D.chem_components)
        eout["num_helical_segments"] = len(hseg)
        eout["num_single-stranded_segments"] = len(sseg)
        eout["helical_segments"] = hseg
        eout["single-stranded_segments"] = sseg
        eout["type"] = etype
        eout["num_nucleotides"] = ntotal
        
        # Update chemical modifications
        updateChemicalModifications(eout["nucleotides"], nuc_dict, eout["chemical_modifications"])
        
        # Get order of strands in each helix
        strand_segments = []
        for h in hseg:
            strand_ids1 = []
            for nid in h["ids1"]:
                sid = strand_dict[nid]["strand_id"]
                if(sid not in strand_ids1):
                    strand_ids1.append(sid)
            strand_ids2 = []
            for nid in reversed(h["ids2"]):
                sid = strand_dict[nid]["strand_id"]
                if(sid not in strand_ids2):
                    strand_ids2.append(sid)
            strand_segments.append(strand_ids1)
            strand_segments.append(strand_ids2)
        
        # Create list of ordered nucleotide ids
        ordered_ids = [] # order of nucleotide IDs which corresponds to DBN
        seen = {}
        for seg in strand_segments:
            for sid in seg:
                if sid in seen:
                    continue
                ordered_ids += strand_dict[sid]["ids"]
                seen[sid] = True
        for s in eout["strands"]:
            if s["strand_id"] not in seen:
                ordered_ids += s["ids"]
        
        if set(ordered_ids) != set(eout["nucleotides"]):
            print(set(eout["nucleotides"]) - set(ordered_ids))
            print("DBN string length and number of nucleotides do not match.", prefix)
        
        # Sort Strands
        si = []
        for s in eout["strands"]:
            m1 = ordered_ids.index(s["3p_end"])
            si.append(m1)
        eout["strands"] = [s for _, s in sorted(zip(si, eout["strands"]))]
        eout["id"] = [c for _, c in sorted(zip(si, eout["id"]))]
        eout["id"] = getHash(*eout["id"])
        out["entities"].append(eout)
        
        # Clean-up files
        #os.remove("{}-dssr.json".format(ename))
        os.remove(ename)
        ecount += 1
    out["num_entities"] = len(out["entities"])
    OUT_DATA.append(out)
    
    return OUT_DATA
