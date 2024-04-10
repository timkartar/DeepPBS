#!/usr/bin/env python
""" This script will serve as main processing script to construct the 
protein/DNA graphs and corresponding features. Output will be a numpy 
.npz file"""

import argparse
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("data_file", help="list of co-crystal structures and corresponding PWM files")
arg_parser.add_argument("config_file", help="list of co-crystal structures and corresponding PWM files")
arg_parser.add_argument("--no_pwm",  action='store_true', default=False)
arg_parser.add_argument("--no_cleanp",  action='store_true', default=False)
ARGS = arg_parser.parse_args()

# built-in modules
import sys
import json
import os
import numpy as np
import random

## should be returned by processDNA
intra_bp_parameters = ["buckle", "shear", "stretch", "stagger", "propeller", "opening"]
inter_bp_parameters = ["shift", "slide", "rise", "tilt", "roll", "twist"]

backbone_parameters = ["major_groove_3dna", "minor_groove_3dna"]
#backbone_parameters = ["major_groove_curves", "minor_groove_curves"]

# deeppbs modules
from deeppbs import StructureData, splitEntities, cleanProtein, processDNA, cleanDNA
from deeppbs import makeDNACG, makeProteinGraph, loadPWM, alignPWMSeq, computeYAndMask
from deeppbs import getAtomSASA, getAtomDepth, getAchtleyFactors, getCV, countContacts

C = json.load(open(ARGS.config_file,"r")) 
outdir = C.get("FEATURE_DATA_PATH", "./output")

for line in [l.strip() for l in open(ARGS.data_file,"r").readlines()]:
    try:
        if line[0] == "#":
            continue
        pdb_file = line.split(",")[0] # this should include file extension. Could be .ent, .pdb, .cif etc...

        if not ARGS.no_pwm:
            pwm_id = line.split(",")[1]
            pwm = loadPWM(pwm_id) # example MA1563.1.jaspar, loaded from ../deeppbs/_data/pwms.pickle

        try:
            structure = StructureData(os.path.join(C["PDB_FILES_PATH"], pdb_file), name="co_crystal")
        except:
            continue


        # split complex into protein and DNA
        protein, dna = splitEntities(structure)

        # clean protein and DNA to fix missing atoms etc...
        
        if not ARGS.no_cleanp:
            try:
                protein, _  = cleanProtein(protein, add_charge_radius=True)
            except Exception as e:
                print("ERROR: clean protein error", pdb_file, e)
            
        dna = cleanDNA(dna,  fix_modified_nucleotide_hetflags=True)
        
        ### DNA-specific features
        # generate DNA CG-bead point cloud representing atom groups
        #try:
        dna_data = processDNA(dna, quiet=False)
        #except Exception as e:
        #    print("ERROR: modified nucleotide", pdb_file, e)
        #    continue

        dna_helices = []
        for entity in dna_data[0]["entities"]:
            dna_helices += entity["helical_segments"]
        try:
            assert len(dna_helices) == 1, "number of helices present is not 1 ! Should have been taken care of in a previous steps"
        except:
            print("ERROR: helix count problem", len(dna_helices), pdb_file)
            continue
        dna_helix = dna_helices[0]


        try:
            V_dna, dna_seq, fn, dna_vectors = makeDNACG(dna, dna_helix) # dna atom group positions
        except Exception as e:
            print("ERROR: missing C5/OP1/OP1/OP2")
            continue

        if not ARGS.no_pwm:
            Y_pwm, pwm_mask, dna_mask, aln_score = computeYAndMask(pwm, dna_seq)
        else:
            Y_pwm = dna_seq
            pwm_mask = [True]*dna_seq.shape[0]
            dna_mask = [True]*dna_seq.shape[0]
            aln_score = [None]
        #print("SCORE:", aln_score, pdb_file)
        
        # add shape information as featues
        N = dna_helix["length"]
        X_dna = np.zeros((N, 14))
        dna_feature_names = []
        i = 0
        # intra-bp shape parameters
        for param in intra_bp_parameters:
            X_dna[:,i] = np.array(dna_helix['shape_parameters'][param])
            i += 1
            dna_feature_names.append(param)

        # inter-bp shape parameters
        for param in inter_bp_parameters:
            s = dna_helix['shape_parameters'][param]
            s.append(np.mean(dna_helix['shape_parameters'][param])) # improved padding values here
            X_dna[:,i] = np.array(s)
            i += 1
            dna_feature_names.append(param)
            '''
            s = dna_helix['shape_parameters'][param]
            s_di = s.copy()
            s = []
            s.append(s_di[0])
            for idx in range(len(s_di)-1):
                 s.append((s_di[idx] + s_di[idx+1])/2)
            s.append(s_di[-1])
            s = np.array(s)
            #print(s, "1", len(s))
            #s.append(0) # TODO: improved padding values here
            X_dna[:,i] = np.array(s)
            i += 1
            dna_feature_names.append(param)
            '''

        # backbone shape parameters
        for param in backbone_parameters:
            s = list(map(lambda x: 0 if x == 'NA' else x, dna_helix['shape_parameters'][param]))
            if '3dna' in param:
                s_di = s.copy()
                s = []
                s.append(s_di[0])
                for idx in range(len(s_di)-1):
                    s.append((s_di[idx] + s_di[idx+1])/2)
                           
                s.append(s_di[-1])
                s = np.array(s)
            X_dna[:,i] = np.array(s)
            i += 1
            dna_feature_names.append(param)
        
        #try:
        #    V_dna.shape
        #except:
        #    print("ERROR: missing atom problem", pdb_file)
        #    continue
        X_dna_point = np.zeros((V_dna.shape[0]*V_dna.shape[1], V_dna.shape[1] + fn.shape[2] + X_dna.shape[1]))
        
        for i in range(V_dna.shape[0]):
            for j in range(V_dna.shape[1]):
                X_dna_point[i*V_dna.shape[1] + j,j] = 1
                X_dna_point[i*V_dna.shape[1] + j, V_dna.shape[1]:V_dna.shape[1] + fn.shape[2]] = fn[i,j,:]
                X_dna_point[i*V_dna.shape[1] + j, V_dna.shape[1]+fn.shape[2]:] = X_dna[i,:]
        
        ### Protein-specific features


        pro_features = ["charge","radius"] # list of feature names
        pro_features.append(getAtomSASA(protein, classifier=None))
        pro_features += getAchtleyFactors(protein)
        pro_features.append(getCV(protein, 7.5, feature_name="cv", impute_hydrogens=True))
        
        V_prot, X_prot, E_prot, prot_vectors =  makeProteinGraph(protein, feature_names=pro_features)
        
        contacts = countContacts(protein, pdb_file, V_dna, dna_mask[0])
        print("CONTACT COUNT", contacts[0], pdb_file)

        if not ARGS.no_pwm:
            np.savez_compressed(os.path.join(outdir,"{}_{}.npz".format(pdb_file.split("/")[-1].replace(".pdb","").replace(".cif",""), pwm_id)), V_dna=V_dna, X_dna=X_dna, X_dna_point=X_dna_point, dna_feature_names=dna_feature_names, V_prot=V_prot, X_prot=X_prot, E_prot=E_prot, prot_feature_names=pro_features, Y_hard=dna_seq,Y_pwm=Y_pwm, pwm_mask=pwm_mask, dna_mask=dna_mask, aln_score=np.array([aln_score]), dna_vectors=dna_vectors, prot_vectors=prot_vectors, contacts=contacts)
        else:
            np.savez_compressed(os.path.join(outdir,"{}.npz".format(pdb_file.split("/")[-1].replace(".pdb","").replace(".cif",""))), V_dna=V_dna, X_dna=X_dna, X_dna_point=X_dna_point, dna_feature_names=dna_feature_names, V_prot=V_prot, X_prot=X_prot, E_prot=E_prot, prot_feature_names=pro_features, Y_hard=dna_seq,Y_pwm=Y_pwm, pwm_mask=pwm_mask, dna_mask=dna_mask, aln_score=np.array([aln_score]), dna_vectors=dna_vectors, prot_vectors=prot_vectors, contacts=contacts)
    except:
        continue
