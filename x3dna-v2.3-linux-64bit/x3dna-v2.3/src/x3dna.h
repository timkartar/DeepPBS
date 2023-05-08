#ifndef _X3DNA_H
#define _X3DNA_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#define NR_END 1  /* for NRC related functions */
#define FREE_ARG char*

#define NO_MATCH -1L
#define DUMMY -1L
#define TRUE 1L
#define FALSE 0L

#define BUF32 32
#define BUF512 512
#define BUF2K 2048
#define BUFBIG 8192

#define UNUSED_PARAMETER(x) (void)(x)
#define NELEMS(x) ((sizeof (x))/(sizeof ((x)[0])))

/* ********* SLRE: http://code.google.com/p/slre/ */
enum slre_option { SLRE_CASE_SENSITIVE = 0, SLRE_CASE_INSENSITIVE = 1 };
enum slre_capture { SLRE_STRING, SLRE_INT, SLRE_FLOAT };

typedef struct {
    long min_base_hb;
    double hb_lower;
    double hb_dist1;
    double hb_dist2;
    char hb_atoms[BUF512];
    long hb_idx[BUF512];
    char alt_list[BUF512];

    double max_dorg;
    double min_dorg;
    double max_dv;
    double min_dv;
    double max_plane_angle;
    double min_plane_angle;
    double max_dNN;
    double min_dNN;

    double helix_break;
    double std_curved;

    double water_dist;
    double water_dlow;
    char water_atoms[BUF512];
    long water_idx[BUF512];

    double o3p_dist;
} miscPars;

typedef struct {
    long DEBUG;
    long VERBOSE;
    long NUM_ELE;
    long CHAIN_CASE;
    long ALL_MODEL;
    long ATTACH_RESIDUE;  /* following Pascal's request */
    long THREE_LETTER_NTS;  /* output 3-letter nucleotides: ADE/CYT/GUA/THY/URA */
    long PDBV3;  /* PDB v3.x; OP1/OP2/C7; DA/DC/DG/DT etc */
    long ORIGINAL_COORDINATE;  /* use original PDB coordinates in bestpairs.pdb etc */
    long OCCUPANCY;  /* if to check for atom occupancy property */
    long HEADER;  /* how much header info to output */

    double NT_CUTOFF;  /* RMSD cutoff for identify a nucleotide */

    char X3DNA_VER[BUF512];
    char X3DNA_HOMEDIR[BUF512];
    char CHAIN_MARKERS[BUF512];
    char REBUILD_CHAIN_IDS[BUF512];

    char *PROGNAME;
    char **ATOM_NAMES;

    long NUM_SATOM;
    char **ATOMLIST;

    long NUM_SBASE;
    char **BASELIST;

    char **AtomName0;
    char **ResName0;
    long Name0;

    miscPars misc_pars;
} struct_Gvars;

/* global variables declaration */
extern struct_Gvars Gvars;

#define DEBUG_LEVEL 6

#define SEPC '\t'  /* tab-delimited fields */
#define WSPACES " ,\t\n"  /* space, comma, tab & newline */
#define SKIPS "#\0"  /* lines to be skipped */

#define UNKATM "XX"
#define DUMSTR "XXXXXX"
#define SFACTOR    2L  /* scale factor for increasing memory */
#define PI 3.141592653589793
#define XEPS 1.0e-7
#define XBIG 1.0e+18
#define XBIG_CUTOFF 1.0e+16
#define MFACTOR 10000.0
#define NMISC 34  /* number of characters of miscellaneous items */
#define NATOMCOL 11  /* number of atom with defined colors */
#define NBASECOL 7  /* number of base with defined colors */

#define PS_DFTSIZE 500  /* 500 points PS default size */
#define PS_BOUND 10  /* boundary offset */
#define FIG_DFTSIZE 8333  /* 8333 units XFIG default size */
#define FIG_BOUND 166  /* boundary offset: 10/72*1200 */

#define PAR_FILE "misc_3dna.par"  /* miscellaneous parameters */
#define BASE_FILE "baselist.dat"  /* 3-letter to 1-letter base residue */
#define ATOM_FILE "atomlist.dat"  /* 2-letter to atomic symbol */
#define HELP3DNA "help3dna.dat"  /* help file for 3DNA */

#define REF_FILE "ref_frames.dat"  /* reference frames */
#define MREF_FILE "mref_frames.dat"  /* multiplet reference frames in <find_pair> */
#define POC_FILE "poc_haxis.r3d"  /* P, O4' and C1' atom radius in r3d format */
#define MUL_FILE "multiplets.pdb"  /* multiplets */
#define ALLP_FILE "allpairs.pdb"  /* all base-pairs */
#define BESTP_FILE "bestpairs.pdb"  /* best base-pairs */
#define STACK_FILE "stacking.pdb"  /* for stacking diagram */
#define HSTACK_FILE "hstacking.pdb"  /* for stacking w.r.t. to middle helical frame */
#define HLXREG_FILE "hel_regions.pdb"  /* helical regions in <find_pair> */
#define BPORDER_FILE "bp_order.dat"  /* base-pair ordering in <find_pair> */
#define COLCHN_FILE "col_chains.scr"  /* color chains in <find_pair> */
#define COLHLX_FILE "col_helices.scr"  /* color helices in <find_pair> */
#define AUX_FILE "auxiliary.par"  /* auxiliary parameters in <analyze> */
#define BPSTEP_FILE "bp_step.par"  /* base-pair step parameters in <analyze> */
#define HLXSTEP_FILE "bp_helical.par"  /* helical step parameters in <analyze> */
#define SEVEN_FILE "cf_7methods.par"  /* compare seven methods in <analyze> */
#define HB_FILE "hbonds_info.dat"  /* H-bonding information <r3d_atom/stack2img> */
#define LKG_FILE "bonds_lkg.dat"  /* bond linkage, as in HB_FILE */
#define SNUM_FILE "serial_num.pdb"  /* PDB file with serial atom numbers */
#define ATOMALC_FILE "atom_lkg.alc"  /* atomic linkage <r3d_atom/stack2img> */
#define BBLKALC_FILE "bblk_lkg.alc"  /* base block linkage <pdb2img> */
#define TMP_FILE "tmp_file"  /* temporary multiplets input file <find_pair> */
#define MULBP_FILE "mulbp.inp"  /* multiplets input file <find_pair> */
#define ROTMAT_FILE "rotmat.dat"  /* rotation matrix <rotate_mol> */
#define VIEW1_FILE "pmiview1"  /* PMI view #1 file <rotate_mol> */
#define VIEW2_FILE "pmiview2"  /* PMI view #2 file <rotate_mol> */
#define VIEW3_FILE "pmiview3"  /* PMI view #3 file <rotate_mol> */

#define NP 101L  /* maximum number of pairs per base */
#define BOND_UPPER_LIMIT 2.5  /* for function torsion */
#define HTWIST0 0.05  /* minimum helical twist */
#define BOND_FACTOR 1.15  /* bond distance criterion */
#define NBOND_FNUM  2.0  /* estimated # of bond from # of atoms */
#define NON_WC_IDX  6  /* non-Watson-Crick base index */
#define AXIS_LENGTH 3.5  /* reference axis length */
#define NUM_BASE_ATOMS BUF512  /* max. no. of base atoms in a residue */
#define NUM_RESIDUE_ATOMS BUF512  /* max. no. of atoms in a residue */
#define NUM_DINUCLEOTIDE_ATOMS BUFBIG  /* max. no. of atoms per dinucleotide */
#define EMPTY_NUMBER -9999.99
#define EMPTY_CRITERION -9999
#define MAXBASE 30000  /* maximum number of bases in regular/fiber */
#define NELE 12  /* 12 elements */
#define O3P_UPPER 2.5  /* upper limit for O3-P connection */
#define RTNNUM 37  /* number of returned values from check_pair */
#define PSTNUM 29  /* number of parameters kept in pair_stat */
#define OLCRT 1.2  /* criterion for overlaps */
#define MBASES 50  /* maximum bases per unit */
#define MAXCH 100  /* maximum # of chains in a biological unit */
#define END_STACK_XANG 125.0  /* find_pair; bdl070 (105 deg) */
#define MAXCLEN 52  /* Maximum length of a color name */

#define RA_LIST " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
#define WC_LIST "XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"
#define CX_LIST "ACGITUX"
#define CB_LIST "ACGITU"
#define NT_LIST "  A", "  C", "  G", "  I", "  T", "  U", \
                "ADE", "CYT", "GUA", "INO", "THY", "URA", \
                " +A", " +C", " +G", " +I", " +T", " +U"

#define WATER_LIST "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP"

#define NUM_SAA 20  /* number of standard amino acids */
#define NUM_ATM 12  /* number of side chain atoms, excluding Ca */
#define AA_LIST "ALA", "VAL", "PHE", "PRO", "MET", "ILE", "LEU",  /* hydrophobic */ \
                "ASP", "GLU", "LYS", "ARG",  /* charged */ \
                "SER", "THR", "TYR", "HIS", "CYS", "ASN", "GLN", "TRP",  /* polar */ \
                "GLY"  /* smallest */

#define OVERLAP 0.01

#define PBLKALC_FILE "pblk_lkg.alc"  /* peptide block linkage */
#define SNAP_PEP_PDB "snap_pep.pdb"
#define SNAP_PEP_ALC "snap_pep.alc"
#define TRSP_RMS 0.25
#define DUMCHAR '@'

/* personal_begin */
/* for snap */
#define LBIG 1000000
#define DNA_BASE "ACGT"
#define SNAP_AAA "snap_aaa.pdb"
#define ABLKALC_FILE "ablk_lkg.alc"  /* amino acid planar moiety block linkage */
#define SNAP_NTS "snap_nts.pdb"
#define SNAP_OPTS "snap_options"
#define WITH_BASE 1  /* base or side-chain */
#define WITH_BKBN 2  /* backbone */
#define WITH_BOTH 3

#define PDBX "PDBx:"
#define O2_STACK "o2_stack.dat"
/* personal_end */

#include "x3dna_fncs.h"

#endif  /* _X3DNA_H */
