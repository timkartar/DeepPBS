
#### Biopython Disordered Atom Fix ####
try:
    import Bio.PDB
    copy = Bio.PDB.Atom.copy
    def myCopy(self):
        shallow = copy.copy(self)
        for child in self.child_dict.values():
            shallow.disordered_add(child.copy())
        return shallow
    Bio.PDB.Atom.DisorderedAtom.copy=myCopy
except ModuleNotFoundError:
    # BioPython isn't found - futher attempts to import will raise an exception
    pass

# built in modules
import logging
import os
import re
import subprocess
from random import choice
from string import ascii_letters

# third party modules
import numpy as np

# geobind modules
from .strip_hydrogens import stripHydrogens
from .data import data
from .structure_data import StructureData
from .run_pdb2pqr import runPDB2PQR

class ResidueMutator(object):    
    def __init__(self, tripeptides=None, components=None, standard_residues=None):
        """ The mutator object takes a non-standard residue or incomplete residue and modifies it
        """
        try:
            from Bio.PDB import PDBParser
            from Bio.SVDSuperimposer import SVDSuperimposer
        except ModuleNotFoundError:
            raise ModuleNotFoundError("BioPython is required for this functionality")
        
        # get defaults if not provided
        if standard_residues is None:
            standard_residues = data.standard_residues
        if tripeptides is None:
            tripeptides = data.tripeptides
        if components is None:
            components = data.chem_components
        self.components = components
        self.candidates = {}
        self.standard_residues = standard_residues
        self.imposer = SVDSuperimposer()
        self.parser = PDBParser(PERMISSIVE=1,QUIET=True)
        
        # build up candidate structures
        for fn in tripeptides:
            structure = self.parser.get_structure("", fn)
            resn = structure[0][" "][2].get_resname()
            self.candidates[resn] = []
            for model in structure:
                self.candidates[resn].append(model[" "][2])
    
    def mutate(self, residue, replace_backbone=True):
        resn = residue.get_resname()
        
        if self.standard(resn):
            # the residue is already a standard residue, here for repair
            parn = resn
        else:
            parn = self.components[resn]['_chem_comp.mon_nstd_parent_comp_id']
            if not self.standard(parn):
                # the parent residue is a nonstandard residue, can't mutate
                return False
        
        if parn not in self.candidates:
            # parent not in candidate structures
            return False
        
        sc_fixed = set(self.components[resn]['side_chain_atoms']) # side chain atoms of fixed residue
        sc_movin = set(self.components[parn]['side_chain_atoms']) # side chain atoms of standard parent
        atom_names = sc_fixed.intersection(sc_movin)
        
        # get list of side chain atoms present in residue
        atom_list = []
        for atom in atom_names:
            if atom in residue:
                atom_list.append(atom)
        
        if len(atom_list) == 0:
            return False
        
        # get side chain atom coordinates
        fixed_coord = np.zeros((len(atom_list), 3))
        for i in range(len(atom_list)):
            fixed_coord[i] = residue[atom_list[i]].get_coord()
        
        # loop over candidates, finding best RMSD
        moved_coord = np.zeros((len(atom_list), 3))
        min_rms = 99999
        rotm = None
        tran = None
        min_candidate = None
        for candidate in self.candidates[parn]:
            for j in range(len(atom_list)):
                moved_coord[j] = candidate[atom_list[j]].get_coord()
            # perfom SVD fitting
            self.imposer.set(fixed_coord, moved_coord)
            self.imposer.run()
            if self.imposer.get_rms() < min_rms:
                min_rms = self.imposer.get_rms()
                rotm, tran = self.imposer.get_rotran()
                min_candidate = candidate
        
        # copy the candidate to a new object
        candidate = min_candidate.copy()
        candidate.transform(rotm, tran)
        stripHydrogens(candidate)
        
        if replace_backbone:
            # replace backbone atoms of candidate
            backbone_atoms = self.components[resn]['main_chain_atoms']
            for atom in backbone_atoms:
                if atom not in residue:
                    continue
                if atom not in candidate:
                    candidate.add(residue[atom].copy())
                candidate[atom].set_coord(residue[atom].get_coord())
        
        return candidate
    
    def standard(self, resname):
        return resname in self.standard_residues
    
    def modified(self, resname):
        if resname in self.standard_residues:
            # it's standard, not modified
            return False
        
        if resname in self.components and '_chem_comp.mon_nstd_parent_comp_id' in self.components[resname]:
            return (
                (resname not in self.standard_residues)
                and
                (self.components[resname]['_chem_comp.mon_nstd_parent_comp_id'] in self.standard_residues)
            )
        else:
            # has no standard parent field - can't be modified
            return False

def heavyAtomCount(residue):
    count = 0
    for atom in residue:
        if atom.element != "H":
            count += 1
    return count

def tempFileName(prefix, ext):
    return "%s.%s" % (prefix + ''.join(choice(ascii_letters) for i in range(15)), ext)

def cleanProtein(structure,
        mutator=None,
        regexes=None,
        hydrogens=True,
        run_pdb2pqr=True,
        quiet=False,
        remove_numerical_chain_id=False,
        method="geobind",
        **kwargs
    ):
    """ Perform any operations needed to modify the structure or sequence of a protein
    chain.
    """
    prefix = structure.name # used for file names
    
    if remove_numerical_chain_id:
        # APBS and TABI-PB does not process numerical chain IDs correctly. This is a work-around
        available_ids = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
        
        # find current chain ids
        taken_ids = set()
        for chain in structure.get_chains():
            cid = chain.get_id()
            taken_ids.add(cid)
        
        # iterate over chains and update
        chain_map = {}
        for chain in structure.get_chains():
            cid = chain.get_id()
            if cid.isnumeric():
                # we want to replace this chain id
                while len(available_ids) > 0:
                    new_id = available_ids.pop()
                    if new_id in taken_ids:
                        continue
                    else:
                        break
                chain_map[cid] = new_id
                chain.id = new_id
            else:
                chain_map[cid] = cid
    
    if method == "geobind":
        # set up needed objects
        if regexes is None:
            regexes = data.regexes
        if mutator is None:
            mutator = ResidueMutator(data.tripeptides, data.chem_components)
        
        # remove non-standard residues
        for chain in structure.get_chains():
            replace = []
            remove = []
            for residue in chain:
                resn = residue.get_resname().strip()
                resid = residue.get_id()
                if resn in data.chem_components and heavyAtomCount(residue)/(data.chem_components[resn]['heavy_atom_count']-1) < 0.6:
                    # too many missing atoms - replace residue
                    replace.append(resid)
                elif mutator.standard(resn):
                    if resid[0] == ' ':
                        continue
                    else:
                        remove.append((resid, "removed HETATM standard residue: %s"))
                elif resn == 'HOH' or resn == 'WAT':
                    remove.append( (resid, None) )
                elif regexes["SOLVENT_COMPONENTS"].search(resn):
                    continue
                elif mutator.modified(resn):
                    replace.append(resid)
                else:
                    remove.append( (resid, "removed unrecognized residue: %s") )
            
            for rid, reason in remove:
                if reason is not None and not quiet:
                    logging.info(reason, chain[rid].get_resname())
                chain.detach_child(rid)
            
            for rid in replace:
                replacement = mutator.mutate(chain[rid])
                if replacement:
                    if not quiet:
                        logging.info("replacing residue %s with %s", chain[rid].get_resname(), replacement.get_resname())
                    replacement.id = rid
                    idx = chain.child_list.index(chain[rid])
                    chain.child_list[idx] = replacement
                else:
                    if not quiet:
                        logging.info("could not perform replacement on %s, removing", chain[rid].get_resname())
                    chain.detach_child(rid)
    elif method == "pdbfixer":
        try:
            from pdbfixer import PDBFixer
            from openmm.app import PDBFile
        except ModuleNotFoundError:
            raise ModuleNotFoundError("The dependencies 'pdbfixer' and 'openmm' are required with option 'method=\"pdbfixer\"'")
        
        # create a temp file
        tmpFile1 = tempFileName(prefix, 'pdb')
        structure.save(tmpFile1)
        
        # run pdbfixer
        fixer = PDBFixer(filename=tmpFile1)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        
        tmpFile2 = tempFileName(prefix, 'pdb')
        PDBFile.writeFile(fixer.topology, fixer.positions, open(tmpFile2, 'w'), keepIds=True)
        
        # load new fixed structure
        structure = StructureData(tmpFile2, name=prefix)
        
        # clean up
        os.remove(tmpFile1)
        os.remove(tmpFile2)
    
    # run PDB2PQR if requested
    if run_pdb2pqr:
        structure, pqrFile = runPDB2PQR(structure, **kwargs)
    
    # remove hydrogens if requested
    if not hydrogens:
        stripHydrogens(structure)
    
    # decide what to return
    rargs = [structure]
    if run_pdb2pqr:
        rargs.append(pqrFile)
    if remove_numerical_chain_id:
        rargs.append(chain_map)
    
    return tuple(rargs)
