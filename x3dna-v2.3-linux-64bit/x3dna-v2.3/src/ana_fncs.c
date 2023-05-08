#include "x3dna.h"

enum { TOR_XXX = 0, TOR_ALPHA, TOR_BETA, TOR_GAMMA, TOR_DELTA, TOR_EPSILON, TOR_ZETA, TOR_CHI,
    TOR_ETA, TOR_THETA, TOR_ETA1, TOR_THETA1, TOR_ETA2, TOR_THETA2,
    TOR_V0, TOR_V1, TOR_V2, TOR_V3, TOR_V4, TOR_TM, TOR_PHASE
};

enum { NT_NUM = 0, NT_P, NT_O5, NT_C5, NT_C4, NT_C3, NT_O3, NT_C2, NT_C1, NT_O4, NT_N, NT_C,
    O3P_LKG
};

char *nt_atoms[] = { "XXXX", " P  ", " O5'", " C5'", " C4'", " C3'", " O3'", " C2'",
    " C1'", " O4'", "N9N1", "C4C2", "OPLK"
};

/* for calculating overlap area between polygons */
#define MNPOLY 1000  /* maximum size of a polygon */

typedef struct {
    double x;
    double y;
} point;

typedef struct {
    point ip;
    point rx;
    point ry;
    long inside;
} vertex;

static double get_chi360(double chi)
{
    return (chi > 0) ? chi : chi + 360;
}

static long in_trans(double chi360)
{
    return dval_in_range(chi360, 165, 315);
}

static void check_chi(double chi, char *fmt, char *bstr, char *ostr)
{
    parcat(ostr, chi, fmt, bstr);

    /* see http://www.fli-leibniz.de/ImgLibDoc/nana/chi.gif */
    if (chi > EMPTY_CRITERION) {
        double chi360;

        chi360 = get_chi360(chi);

        if (dval_in_range(chi360, 45, 95))  /* [60..80], with 15 allowance */
            strcat(ostr, " syn  ");
        else if (in_trans(chi360))  /* [-60, 180], with 15 allowance: [165, 315] */
            strcat(ostr, " anti ");
        else
            strcat(ostr, "      ");
    } else
        strcat(ostr, " ---  ");
}

/* http://www.imb-jena.de/Piet/help/backbone.html
 *   BI: (epsilon-zeta)	= -160 ... +20
 *  BII: (epsilon-zeta)	=  +20 ... +200 */
static void check_BI_BII(double epsilon, double zeta, char *bstr, char *fmt, char *ostr)
{
    if (epsilon < EMPTY_CRITERION || zeta < EMPTY_CRITERION) {
        parcat(ostr, EMPTY_NUMBER, fmt, bstr);
        strcat(ostr, "  ---");  /* no data for classification */

    } else {
        double d;

        if (epsilon < 0)
            epsilon += 360;
        if (zeta < 0)
            zeta += 360;

        d = epsilon - zeta;
        parcat(ostr, d, fmt, bstr);

        if (d < 0)
            d += 360;

        strcat(ostr, dval_in_range(d, 20, 200) ? "  BII" : "  BI");
    }
}

static void output_header_sugar_torsion(FILE * fp)
{
    fprintf(fp, "Sugar conformational parameters: \n\n"
            "Note: v0: C4'-O4'-C1'-C2'\n"
            "      v1: O4'-C1'-C2'-C3'\n"
            "      v2: C1'-C2'-C3'-C4'\n"
            "      v3: C2'-C3'-C4'-O4'\n"
            "      v4: C3'-C4'-O4'-C1'\n\n"
            "      tm: the amplitude of pucker\n"
            "      P:  the phase angle of pseudorotation\n\n");
}

static void output_header_ss_Zp_Dp(FILE * fp)
{
    fprintf(fp,
            "    ssZp: single-stranded Zp, defined as the z-coordinate of the 3'\n"
            "            phosphorus atom (P) expressed in the standard reference\n"
            "            frame of preceding base; the value is POSITIVE when P lies\n"
            "            on the +z-axis side (base in anti conformation); NEGATIVE\n"
            "            if P is on the -z-axis side (base in syn conformation)\n"
            "      Dp: perpendicular distance of the 3' P atom to the glycosydic bond\n"
            "            [as per the MolProbity paper of Richardson et al. (2010)]\n\n");
}

static void add_sugar_pucker(double phase_angle, char *bstr, char *ostr)
{
    static char *sugar_pucker[10] = { "C3'-endo", "C4'-exo ", "O4'-endo", "C1'-exo ", "C2'-endo",
        "C3'-exo ", "C4'-endo", "O4'-exo ", "C1'-endo", "C2'-exo "
    };
    char temp[BUF512];
    long m;

    if (phase_angle > EMPTY_CRITERION) {
        m = (long) floor(phase_angle / 36.0);
        sprintf(temp, "%12s", sugar_pucker[m]);
        strcat(ostr, temp);
    } else
        strcat(ostr, bstr);
}

static void output_header_bb_torsion(FILE * fp)
{
    fprintf(fp, "Main chain and chi torsion angles: \n\n"
            "Note: alpha:   O3'(i-1)-P-O5'-C5'\n"
            "      beta:    P-O5'-C5'-C4'\n"
            "      gamma:   O5'-C5'-C4'-C3'\n"
            "      delta:   C5'-C4'-C3'-O3'\n"
            "      epsilon: C4'-C3'-O3'-P(i+1)\n"
            "      zeta:    C3'-O3'-P(i+1)-O5'(i+1)\n\n"
            "      chi for pyrimidines(Y): O4'-C1'-N1-C2\n"
            "          chi for purines(R): O4'-C1'-N9-C4\n\n");
}

static void output_header_BI_BII_chi_syn_anti(FILE * fp)
{
    fprintf(fp,
            "          chi in [165, -45(315)] for anti conformation\n"
            "                 chi in [45, 95] for syn conformation\n\n"
            "          e-z: epsilon - zeta\n"
            "              BI:  e-z = [-160, +20]\n" "              BII: e-z = [+20, +200]\n\n");
}

static double get_ss_Zp(double *P_xyz, double *o_xyz, double *n_xyz)
{
    double dd[4];

    ddxyz(o_xyz, P_xyz, dd);
    return dot(dd, n_xyz);
}

static double idx2torsion(long *idx, double **xyz, long chk_lkg)
{
    double d = EMPTY_NUMBER, **xyz4;
    long i;

    xyz4 = dmatrix(1, 4, 1, 3);

    for (i = 1; i <= 4; i++) {
        if (!idx[i])
            break;
        cpxyz(xyz[idx[i]], xyz4[i]);
    }
    if (i > 4)
        d = chk_lkg ? torsion(xyz4) : torsion2(xyz4);

    free_dmatrix(xyz4, 1, 4, 1, 3);

    return d;
}

void populate_nt_info(long num_residue, long **seidx, char **ResName, char *ChainID,
                      long *ResSeq, char **Miscs, char *bseq, char **nt_info)
{
    char idmsg[BUF512];
    long i, ib;

    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        base_str(ChainID[ib], ResSeq[ib], Miscs[ib], ResName[ib], bseq[i], 1, idmsg);
        if (str_pmatch(idmsg, "....>"))
            strcpy(nt_info[i], idmsg + 5);
        else
            strcpy(nt_info[i], idmsg);
    }
}

void populate_nt_list(long num_residue, long **seidx, long *RY, char *bseq,
                      char **AtomName, double **xyz, long **nt_list)
{
    long i, j, k, ib, ie;

    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];

        for (j = 1; j <= NT_O4; j++)
            nt_list[i][j] = find_1st_atom(nt_atoms[j], AtomName, ib, ie, "");

        if (RY[i] == 1) {  /* R: puRine */
            nt_list[i][NT_N] = find_1st_atom(" N9 ", AtomName, ib, ie, "");
            nt_list[i][NT_C] = find_1st_atom(" C4 ", AtomName, ib, ie, "");
        } else if (RY[i] == 0) {  /* Y: pyRimidine */
            if (bseq[i] == 'P' || bseq[i] == 'p') {  /* pseudo-U */
                nt_list[i][NT_N] = find_1st_atom(" C5 ", AtomName, ib, ie, "");
                nt_list[i][NT_C] = find_1st_atom(" C4 ", AtomName, ib, ie, "");
            } else {
                nt_list[i][NT_N] = find_1st_atom(" N1 ", AtomName, ib, ie, "");
                nt_list[i][NT_C] = find_1st_atom(" C2 ", AtomName, ib, ie, "");
            }
        }

        k = 0;
        for (j = 1; j <= NT_C; j++)
            if (nt_list[i][j])
                k++;
        nt_list[i][NT_NUM] = k;
    }

    for (i = 1; i < num_residue; i++) {
        ib = nt_list[i][NT_O3];
        ie = nt_list[i + 1][NT_P];
        if (ib && ie && within_limits(xyz[ib], xyz[ie], 0.8, O3P_UPPER))
            nt_list[i][O3P_LKG] = TRUE;
    }
}

/* Get PDB data file name, output file name, and pairing information
   pair_num matrix is allocated here but de-allocated elsewhere */
long **read_input(char *inpfile, char *pdbfile, char *outfile, long *ds, long *num_bp,
                  long *ip, long *hetatm)
{
    char str[BUF512];
    long i, j = 0, k, **pair_num;
    FILE *fp;

    fp = open_file(inpfile, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%s", pdbfile) != 1)
        fatal("cannot read pdbfile name: %s\n", pdbfile);  /* PDB file name */

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%s", outfile) != 1)
        fatal("cannot read outfile name: %s\n", outfile);  /* output file name */

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", ds) != 1)
        fatal("cannot read strand information\n");  /* ds is a pointer */
    if (*ds != 2 && *ds != 1)
        fatal("allowed options: 2--duplex and 1--single helix\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", num_bp) != 1)
        fatal("cannot read number of base-pairs\n");  /* num_bp is a pointer */
    if (*num_bp <= 0)
        fatal("illegal number (%ld <= 0) of base-pairs\n", *num_bp);

    if (fgets(str, sizeof str, fp) == NULL ||  /* ip/hetatm are pointers */
        ((j = sscanf(str, "%ld %ld", ip, hetatm)) != 1 && j != 2)) {
        if (feof(fp)) {
            *ip = 0;
            *hetatm = 0;
        } else
            fatal("cannot read integer pair numbering/hetero indicator)\n");
    }

    if (j == 1)
        *hetatm = 0;
    pair_num = lmatrix(1, *ds + 1, 1, *num_bp);

    if (*ip) {  /* user input */
        for (i = 1; i <= *num_bp; i++) {
            if (fgets(str, sizeof str, fp) == NULL)
                fatal("cannot read pair list\n");
            if (*ds == 2) {  /* duplex */
                k = sscanf(str, "%ld %ld %ld", &pair_num[1][i], &pair_num[2][i], &pair_num[3][i]);
                if (k != 2 && k != 3)
                    fatal("two serial numbers required\n");
                if (k == 2 || (k == 3 && !pair_num[3][i] && pair_num[3][i] != 1
                               && pair_num[3][i] != 9))
                    pair_num[3][i] = 0;
                if (pair_num[1][i] <= 0 || pair_num[2][i] <= 0)
                    fatal("residue serial number %ld or %ld <= 0\n", pair_num[1][i],
                          pair_num[2][i]);
            } else {  /* single helix */
                if (sscanf(str, "%ld", &pair_num[1][i]) != 1)
                    fatal("one serial number required: %s\n", str);
                if (pair_num[1][i] <= 0)
                    fatal("residue serial number %ld <= 0\n", pair_num[1][i]);
            }
        }
    }

    close_file(fp);

    return pair_num;
}

void print_header(long ds, long num_bp, long num, char *pdbfile, FILE * fp)
{
    time_t run_time;

    fprintf(fp, "    %s\n", Gvars.X3DNA_VER);
    print_sep(fp, '*', 76);
    fprintf(fp, "1. The list of the parameters given below correspond to"
            " the 5' to 3' direction\n   of strand I and 3' to 5' direction" " of strand II.\n\n");
    fprintf(fp, "2. All angular parameters, except for the phase angle"
            " of sugar pseudo-\n   rotation, are measured in degrees in"
            " the range of [-180, +180], and all\n"
            "   displacements are measured in Angstrom units.\n");

    print_sep(fp, '*', 76);
    fprintf(fp, "File name: %s\n", pdbfile);

    run_time = time(NULL);
    fprintf(fp, "Date and time: %s\n", ctime(&run_time));

    fprintf(fp, "Number of %s%s: %ld\n", (ds == 2) ? "base-pair" : "base",
            (num_bp == 1) ? "" : "s", num_bp);
    fprintf(fp, "Number of atoms: %ld\n", num);

    print_sep(fp, '*', 76);
    print_pdb_title(pdbfile, "*", fp);
}

static void output_atom_xyz(FILE * fp, long idx, double **xyz, char *bstr, char *fmt)
{
    long i;

    if (idx)
        for (i = 1; i <= 3; i++)
            fprintf(fp, fmt, xyz[idx][i]);
    else
        for (i = 1; i <= 3; i++)
            fprintf(fp, "%s", bstr);
}

void output_Borg_P_C1_C4(long num_residue, double **org, double **xyz, long **nt_list,
                         char **nt_info)
{
    char *bstr = "    ---- ", *fmt = "%9.3f";
    long i, j;
    FILE *fp;

    fp = open_file("Borg_P_C1_C4.dat", "w");
    print_sep(fp, '*', 76);
    fprintf(fp, "xyz coordinates of the base origin, the P, C1' and C4' atoms\n\n");
    fprintf(fp,
            "              base      Ox       Oy       Oz       Px       Py       Pz"
            "      C1x      C1y      C1z      C4x      C4y      C4z\n");

    for (i = 1; i <= num_residue; i++) {
        if (nt_list[i][NT_NUM] < 3)
            continue;

        fprintf(fp, "%4ld %s", i, nt_info[i]);
        for (j = 1; j <= 3; j++)
            fprintf(fp, fmt, org[i][j]);
        output_atom_xyz(fp, nt_list[i][NT_P], xyz, bstr, fmt);
        output_atom_xyz(fp, nt_list[i][NT_C1], xyz, bstr, fmt);
        output_atom_xyz(fp, nt_list[i][NT_C4], xyz, bstr, fmt);
        fprintf(fp, "\n");
    }

    close_file(fp);
}

/* Indexes for sugar, chi torsion, C6-C8 etc atoms */
void atom_list(long ds, long num_bp, long **pair_num, long **seidx, long *RY,
               char **bp_seq, char **AtomName, char **ResName, char *ChainID,
               long *ResSeq, char **Miscs, long **phos, long **c6_c8, long **sugar, long **chi)
{
    char c2c4[5], c6c8[5], idmsg[BUF512], n1n9[5];
    long i, ib, ie, ioffset, j, c1, o4, rnum;

    for (i = 1; i <= ds; i++) {
        for (j = 1; j <= num_bp; j++) {
            rnum = pair_num[i][j];
            ib = seidx[rnum][1];
            ie = seidx[rnum][2];

            if (RY[rnum] < 0)
                fatal("Non-base residue: %s\n", ResName[ib]);
            get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);

            /* backbone: P-O5'-C5'-C4'-C3'-O3' */
            phos[i][j] = find_1st_atom(" P  ", AtomName, ib, ie, idmsg);
            phos[i + 2][j] = find_1st_atom(" O1P", AtomName, ib, ie, idmsg);
            phos[i + 4][j] = find_1st_atom(" O2P", AtomName, ib, ie, idmsg);

            /* sugar: C4'-O4'-C1'-C2'-C3' */
            ioffset = (j - 1) * 5;
            sugar[i][ioffset + 1] = find_1st_atom(" C4'", AtomName, ib, ie, idmsg);
            sugar[i][ioffset + 5] = find_1st_atom(" C3'", AtomName, ib, ie, idmsg);
            o4 = find_1st_atom(" O4'", AtomName, ib, ie, idmsg);
            c1 = find_1st_atom(" C1'", AtomName, ib, ie, idmsg);
            sugar[i][ioffset + 2] = o4;
            sugar[i][ioffset + 3] = c1;
            sugar[i][ioffset + 4] = find_1st_atom(" C2'", AtomName, ib, ie, idmsg);

            /* chi(R): O4'-C1'-N9-C4; chi(Y): O4'-C1'-N1-C2 */
            ioffset = (j - 1) * 4;
            chi[i][ioffset + 1] = o4;
            chi[i][ioffset + 2] = c1;
            if (RY[rnum] == 1) {
                strcpy(n1n9, " N9 ");
                strcpy(c2c4, " C4 ");
                strcpy(c6c8, " C8 ");
            } else if (RY[rnum] == 0) {
                strcpy(n1n9, " N1 ");
                strcpy(c2c4, " C2 ");
                strcpy(c6c8, " C6 ");
                if (bp_seq[i][j] == 'P' || bp_seq[i][j] == 'p') {
                    strcpy(n1n9, " C5 ");
                    strcpy(c2c4, " C4 ");
                }
            }
            chi[i][ioffset + 3] = find_1st_atom(n1n9, AtomName, ib, ie, idmsg);
            chi[i][ioffset + 4] = find_1st_atom(c2c4, AtomName, ib, ie, idmsg);
            c6_c8[i][j] = find_1st_atom(c6c8, AtomName, ib, ie, idmsg);
        }
    }
}

static void get_mc_6_torsions(long num_residue, long *RY, long **bb, double **xyz,
                              double **nt_bb_torsion)
{
    long i, idx[5];

    for (i = 1; i <= num_residue; i++) {
        if (RY[i] < 0)  /* not a nt */
            continue;

        /* alpha: O3'(i - 1) - P(i) - O5'(i) - C5'(i) */
        idx[1] = (i > 1) ? bb[i - 1][6] : 0;
        idx[2] = bb[i][1];
        idx[3] = bb[i][2];
        idx[4] = bb[i][3];
        nt_bb_torsion[i][1] = idx2torsion(idx, xyz, TRUE);

        /* beta: P(i) - O5'(i) - C5'(i) - C4'(i) */
        idx[1] = bb[i][1];
        idx[2] = bb[i][2];
        idx[3] = bb[i][3];
        idx[4] = bb[i][4];
        nt_bb_torsion[i][2] = idx2torsion(idx, xyz, TRUE);

        /* gamma: O5'(i) - C5'(i) - C4'(i) - C3'(i) */
        idx[1] = bb[i][2];
        idx[2] = bb[i][3];
        idx[3] = bb[i][4];
        idx[4] = bb[i][5];
        nt_bb_torsion[i][3] = idx2torsion(idx, xyz, TRUE);

        /* delta: C5'(i) - C4'(i) - C3'(i) - O3'(i) */
        idx[1] = bb[i][3];
        idx[2] = bb[i][4];
        idx[3] = bb[i][5];
        idx[4] = bb[i][6];
        nt_bb_torsion[i][4] = idx2torsion(idx, xyz, TRUE);

        /* epsilon: C4'(i) - C3'(i) - O3'(i) - P(i + 1) */
        idx[1] = bb[i][4];
        idx[2] = bb[i][5];
        idx[3] = bb[i][6];
        idx[4] = (i < num_residue) ? bb[i + 1][1] : 0;
        nt_bb_torsion[i][5] = idx2torsion(idx, xyz, TRUE);

        /* zeta: C3'(i) - O3'(i) - P(i + 1) - O5'(i + 1) */
        idx[1] = bb[i][5];
        idx[2] = bb[i][6];
        idx[3] = (i < num_residue) ? bb[i + 1][1] : 0;
        idx[4] = (i < num_residue) ? bb[i + 1][2] : 0;
        nt_bb_torsion[i][6] = idx2torsion(idx, xyz, TRUE);
    }
}

void get_nt_torsion(long num_residue, double **org, double **xyz, long **nt_list,
                    double **nt_torsion)
{
    static double Pconst;
    double phase_angle, **temp_xyz;
    long i, j, k, im1, ip1, idx[BUF32];

    Pconst = sin(PI / 5) + sin(PI / 2.5);
    temp_xyz = dmatrix(1, BUF32, 1, 3);

    init_dmatrix(nt_torsion, 1, num_residue, 1, BUF32, EMPTY_NUMBER);

    for (i = 1; i <= num_residue; i++) {
        im1 = i - 1;
        ip1 = i + 1;

        /* alpha: O3'(i-1) - P(i) - O5'(i) - C5'(i) */
        idx[1] = (i > 1) ? nt_list[im1][NT_O3] : 0;
        idx[2] = nt_list[i][NT_P];
        idx[3] = nt_list[i][NT_O5];
        idx[4] = nt_list[i][NT_C5];
        idx[5] = nt_list[i][NT_C4];  /* beta: P(i) - O5'(i) - C5'(i) - C4'(i) */
        idx[6] = nt_list[i][NT_C3];  /* gamma: O5'(i) - C5'(i) - C4'(i) - C3'(i) */
        idx[7] = nt_list[i][NT_O3];  /* delta: C5'(i) - C4'(i) - C3'(i) - O3'(i) */
        idx[8] = (i < num_residue) ? nt_list[ip1][NT_P] : 0;  /* epsilon: C4'(i) - C3'(i) - O3'(i) - P(i+1) */
        idx[9] = (i < num_residue) ? nt_list[ip1][NT_O5] : 0;  /* zeta: C3'(i) - O3'(i) - P(i+1) - O5'(i+1) */

        /* no need to check for O3P_LKG: implicit with TRUE option below */
        nt_torsion[i][TOR_ALPHA] = idx2torsion(idx, xyz, TRUE);
        nt_torsion[i][TOR_BETA] = idx2torsion(idx + 1, xyz, TRUE);
        nt_torsion[i][TOR_GAMMA] = idx2torsion(idx + 2, xyz, TRUE);
        nt_torsion[i][TOR_DELTA] = idx2torsion(idx + 3, xyz, TRUE);
        nt_torsion[i][TOR_EPSILON] = idx2torsion(idx + 4, xyz, TRUE);
        nt_torsion[i][TOR_ZETA] = idx2torsion(idx + 5, xyz, TRUE);

        /* chi: O4'-C1'-N9-C4 (R); O4'-C1'-N1-C2 (Y) */
        idx[1] = nt_list[i][NT_O4];
        idx[2] = nt_list[i][NT_C1];
        idx[3] = nt_list[i][NT_N];
        idx[4] = nt_list[i][NT_C];
        nt_torsion[i][TOR_CHI] = idx2torsion(idx, xyz, TRUE);

        /* eta:    C4'(i-1)-P(i)-C4'(i)-P(i+1)    --- C4'
           theta:  P(i)-C4'(i)-P(i+1)-C4'(i+1)
           eta':   C1'(i-1)-P(i)-C1'(i)-P(i+1)    --- C1'
           theta': P(i)-C1'(i)-P(i+1)-C1'(i+1)
           eta":   O(i-1)-P(i)-O(i)-P(i+1)        --- Origin of base
           theta": P(i)-O(i)-P(i+1)-O(i+1)  */
        if (i < num_residue && nt_list[i][O3P_LKG]) {
            idx[1] = (i > 1 && nt_list[im1][O3P_LKG]) ? nt_list[im1][NT_C4] : 0;
            idx[2] = nt_list[i][NT_P];
            idx[3] = nt_list[i][NT_C4];
            idx[4] = nt_list[ip1][NT_P];
            idx[5] = nt_list[ip1][NT_C4];
            nt_torsion[i][TOR_ETA] = idx2torsion(idx, xyz, FALSE);
            nt_torsion[i][TOR_THETA] = idx2torsion(idx + 1, xyz, FALSE);

            /* using C1' instead of C4' */
            idx[1] = (i > 1 && nt_list[im1][O3P_LKG]) ? nt_list[im1][NT_C1] : 0;
            idx[3] = nt_list[i][NT_C1];
            idx[5] = nt_list[ip1][NT_C1];
            nt_torsion[i][TOR_ETA1] = idx2torsion(idx, xyz, FALSE);
            nt_torsion[i][TOR_THETA1] = idx2torsion(idx + 1, xyz, FALSE);

            /* using origin of base reference frame */
            if (idx[2] && idx[4]) {
                if (i > 1)
                    copy_dvector(temp_xyz[1], org[im1], 1, 3);
                copy_dvector(temp_xyz[2], xyz[idx[2]], 1, 3);
                copy_dvector(temp_xyz[3], org[i], 1, 3);
                copy_dvector(temp_xyz[4], xyz[idx[4]], 1, 3);
                copy_dvector(temp_xyz[5], org[ip1], 1, 3);
                if (i > 1)
                    nt_torsion[i][TOR_ETA2] = torsion2(temp_xyz);
                nt_torsion[i][TOR_THETA2] = torsion2(temp_xyz + 1);
            }
        }

        /* sugar ring torsion angles v0 to v4 */
        idx[1] = nt_list[i][NT_C4];
        idx[2] = nt_list[i][NT_O4];
        idx[3] = nt_list[i][NT_C1];
        idx[4] = nt_list[i][NT_C2];
        idx[5] = nt_list[i][NT_C3];
        idx[6] = nt_list[i][NT_C4];
        idx[7] = nt_list[i][NT_O4];
        idx[8] = nt_list[i][NT_C1];
        k = 0;
        for (j = 0; j <= 4; j++) {
            nt_torsion[i][TOR_V0 + j] = idx2torsion(idx + j, xyz, TRUE);
            if (idx[j + 1])
                k++;
        }
        if (k != 5)
            continue;

        /* phase angle and amplitude of pseudorotation */
        phase_angle = atan2(nt_torsion[i][TOR_V4] + nt_torsion[i][TOR_V1] -
                            nt_torsion[i][TOR_V3] - nt_torsion[i][TOR_V0],
                            2.0 * nt_torsion[i][TOR_V2] * Pconst);
        nt_torsion[i][TOR_TM] = nt_torsion[i][TOR_V2] / cos(phase_angle);

        phase_angle = rad2deg(phase_angle);
        if (phase_angle < 0)
            phase_angle += 360;
        nt_torsion[i][TOR_PHASE] = phase_angle;
    }

    free_dmatrix(temp_xyz, 1, BUF32, 1, 3);
}

void get_ss_Zp_Dp(long num_residue, double **org, double **orien, double **xyz,
                  long **nt_list, double **ss_Zp_Dp)
{
    long i, ip1, idxP, idxN, idxC1;

    init_dmatrix(ss_Zp_Dp, 1, 2, 1, num_residue, EMPTY_NUMBER);

    for (i = 1; i <= num_residue; i++) {
        if (nt_list[i][NT_NUM] < 3)
            continue;

        ip1 = i + 1;
        if (i < num_residue && nt_list[i][O3P_LKG] && nt_list[ip1][NT_P]) {
            idxP = nt_list[ip1][NT_P];
            ss_Zp_Dp[1][i] = get_ss_Zp(xyz[idxP], org[i], orien[i] + 6);

            idxN = nt_list[i][NT_N];
            idxC1 = nt_list[i][NT_C1];
            if (idxC1 && idxN)
                ss_Zp_Dp[2][i] = get_point2line_perp_distance(xyz[idxP], xyz[idxC1], xyz[idxN]);
        }
    }
}

static long has_model_number(long num_residue, char **nt_info, long **nt_list)
{
    long i;

    for (i = 1; i <= num_residue; i++) {
        if (nt_list[i][NT_NUM] < 3)
            continue;
        if (nt_info[i][4] == '>')
            return TRUE;
    }

    return FALSE;
}

void output_nt_torsion(long num_residue, char **nt_info, long **nt_list,
                       double **nt_torsion, double **ss_Zp_Dp, FILE * fp)
{
    char *bstr = "    --- ", *fmt = "%8.1f";
    char ostr[BUF512];
    long i, j, with_model_number;

    with_model_number = has_model_number(num_residue, nt_info, nt_list);

    print_sep(fp, '*', 76);
    output_header_bb_torsion(fp);
    output_header_BI_BII_chi_syn_anti(fp);
    fprintf(fp, "              %sbase      chi A/S     alpha    beta   gamma   delta"
            "  epsilon   zeta     e-z BI/BII\n", with_model_number ? "     " : "");
    for (i = 1; i <= num_residue; i++) {
        if (nt_list[i][NT_NUM] < 3)
            continue;

        sprintf(ostr, "%4ld %s", i, nt_info[i]);
        check_chi(nt_torsion[i][TOR_CHI], fmt, bstr, ostr);

        for (j = TOR_ALPHA; j <= TOR_ZETA; j++)
            parcat(ostr, nt_torsion[i][j], fmt, bstr);
        check_BI_BII(nt_torsion[i][TOR_EPSILON], nt_torsion[i][TOR_ZETA], bstr, fmt, ostr);
        fprintf(fp, "%s\n", ostr);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "Pseudo (virtual) eta/theta torsion angles:\n\n");
    fprintf(fp,
            "Note: eta:    C4'(i-1)-P(i)-C4'(i)-P(i+1)\n"
            "      theta:  P(i)-C4'(i)-P(i+1)-C4'(i+1)\n\n"
            "      eta':   C1'(i-1)-P(i)-C1'(i)-P(i+1)\n"
            "      theta': P(i)-C1'(i)-P(i+1)-C1'(i+1)\n\n"
            "      eta\":   Borg(i-1)-P(i)-Borg(i)-P(i+1)\n"
            "      theta\": P(i)-Borg(i)-P(i+1)-Borg(i+1)\n\n");
    fprintf(fp,
            "              %sbase      eta   theta    eta'  theta'    eta\"  theta\"\n",
            with_model_number ? "     " : "");

    for (i = 1; i <= num_residue; i++) {
        if (nt_list[i][NT_NUM] < 3)
            continue;

        sprintf(ostr, "%4ld %s", i, nt_info[i]);
        for (j = TOR_ETA; j <= TOR_THETA2; j++)
            parcat(ostr, nt_torsion[i][j], fmt, bstr);
        fprintf(fp, "%s\n", ostr);
    }

    print_sep(fp, '*', 76);
    output_header_sugar_torsion(fp);
    output_header_ss_Zp_Dp(fp);
    fprintf(fp, "              %sbase       v0      v1      v2      v3      v4"
            "     tm       P    Puckering    ssZp     Dp\n", with_model_number ? "     " : "");
    for (i = 1; i <= num_residue; i++) {
        if (nt_list[i][NT_NUM] < 3)
            continue;

        sprintf(ostr, "%4ld %s", i, nt_info[i]);
        for (j = TOR_V0; j <= TOR_PHASE; j++)
            parcat(ostr, nt_torsion[i][j], fmt, bstr);
        add_sugar_pucker(nt_torsion[i][TOR_PHASE], bstr, ostr);

        if (nt_torsion[i][TOR_PHASE] < EMPTY_CRITERION)
            strcat(ostr, "    ");
        parcat(ostr, ss_Zp_Dp[1][i], "%8.2f", bstr);
        parcat(ostr, ss_Zp_Dp[2][i], "%8.2f", bstr);

        fprintf(fp, "%s\n", ostr);
    }
}

void get_nt_bb_torsion(double **nt_bb_torsion, long num_residue, long **seidx,
                       long *RY, char **AtomName, char **ResName, char *ChainID,
                       long *ResSeq, char **Miscs, double **xyz)
{
    char idmsg[BUF512];
    char *bb_atoms[] = { " P  ", " O5'", " C5'", " C4'", " C3'", " O3'", " C1'" };
    long i, j, ib, ie, num_bb = sizeof(bb_atoms) / sizeof(bb_atoms[0]);
    long **bb;

    bb = lmatrix(1, num_residue, 1, num_bb);
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] < 0)
            continue;

        ib = seidx[i][1];
        ie = seidx[i][2];
        get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);
        for (j = 1; j <= num_bb; j++)
            bb[i][j] = find_1st_atom(bb_atoms[j - 1], AtomName, ib, ie, idmsg);
    }

    init_dmatrix(nt_bb_torsion, 1, num_residue, 1, 6, EMPTY_NUMBER);
    get_mc_6_torsions(num_residue, RY, bb, xyz, nt_bb_torsion);

    free_lmatrix(bb, 1, num_residue, 1, num_bb);
}

static void get_chi_torsions(long ds, long num_bp, long **chi, double **xyz, double **chi_angle)
{
    double **xyz4;
    long i, j, idx, ioffset, m;

    xyz4 = dmatrix(1, 4, 1, 3);

    for (i = 1; i <= ds; i++) {
        for (j = 1; j <= num_bp; j++) {
            ioffset = (j - 1) * 4;
            for (m = 1; m <= 4; m++) {
                idx = chi[i][ioffset + m];
                if (!idx)
                    break;
                cpxyz(xyz[idx], xyz4[m]);
            }
            if (m > 4)  /* all 4 indexes are okay  */
                chi_angle[i][j] = torsion(xyz4);
        }
    }

    free_dmatrix(xyz4, 1, 4, 1, 3);
}

static void get_sugar_torsions(long ds, long num_bp, long **sugar, double **xyz,
                               double **sugar_angle)
{
    static double Pconst;
    static long vidx[5][4] =  /* index of v0 to v4 */
    {
        {1, 2, 3, 4},
        {2, 3, 4, 5},
        {3, 4, 5, 1},
        {4, 5, 1, 2},
        {5, 1, 2, 3}
    };
    double P_angle, **xyz4;
    long i, j, k, i5, idx, ioffset, m, o7;

    Pconst = sin(PI / 5) + sin(PI / 2.5);
    xyz4 = dmatrix(1, 4, 1, 3);

    for (i = 1; i <= ds; i++) {
        for (j = 1; j <= num_bp; j++) {
            /* sugar ring torsion angles v0 to v4 */
            i5 = 0;
            ioffset = (j - 1) * 5;
            o7 = (j - 1) * 7;
            for (k = 1; k <= 5; k++) {
                for (m = 1; m <= 4; m++) {
                    idx = sugar[i][ioffset + vidx[k - 1][m - 1]];
                    if (!idx)
                        break;
                    cpxyz(xyz[idx], xyz4[m]);
                }
                if (m > 4) {  /* all 4 indexes are okay  */
                    sugar_angle[i][o7 + k] = torsion(xyz4);
                    i5++;
                }
            }

            /* phase angle and amplitude of pseudorotation */
            if (i5 == 5) {
                P_angle = atan2(sugar_angle[i][o7 + 5] + sugar_angle[i][o7 + 2]
                                - sugar_angle[i][o7 + 4] - sugar_angle[i][o7 + 1],
                                2 * sugar_angle[i][o7 + 3] * Pconst);
                sugar_angle[i][o7 + 6] = sugar_angle[i][o7 + 3] / cos(P_angle);
                P_angle = rad2deg(P_angle);
                if (P_angle < 0)
                    P_angle = P_angle + 360;
                sugar_angle[i][o7 + 7] = P_angle;
            }
        }
    }

    free_dmatrix(xyz4, 1, 4, 1, 3);
}

static void print_duplex_torsions(long num_bp, long **pair_num, char **bp_seq,
                                  double **chi_angle, double **nt_bb_torsion,
                                  char *bstr, char *fmt, FILE * fp)
{
    char str[BUF512];
    long i, j, m, idx, ds = 2;

    output_header_bb_torsion(fp);

    for (i = 1; i <= ds; i++) {
        (i == 1) ? fprintf(fp, "Strand I\n") : fprintf(fp, "Strand II\n");
        fprintf(fp, "  base    alpha    beta   gamma   delta  epsilon   zeta    chi\n");

        for (j = 1; j <= num_bp; j++) {
            sprintf(str, "%4ld %c ", j, bp_seq[i][j]);
            idx = pair_num[i][j];  /* residue number */
            for (m = 1; m <= 6; m++)
                parcat(str, nt_bb_torsion[idx][m], fmt, bstr);
            parcat(str, chi_angle[i][j], fmt, bstr);
            fprintf(fp, "%s\n", str);
        }

        if (i == 1)
            fprintf(fp, "\n");
    }
}

static void print_ss_torsions(long num_bp, long **pair_num, char **bp_seq,
                              double **chi_angle, double **nt_bb_torsion,
                              char *bstr, char *fmt, FILE * fp)
{
    char str[BUF512];
    long i = 1, j, m, idx;

    output_header_bb_torsion(fp);

    fprintf(fp, "  base    alpha    beta   gamma   delta  epsilon   zeta    chi\n");

    for (j = 1; j <= num_bp; j++) {
        sprintf(str, "%4ld %c ", j, bp_seq[i][j]);
        idx = pair_num[i][j];  /* residue number */
        for (m = 1; m <= 6; m++)
            parcat(str, nt_bb_torsion[idx][m], fmt, bstr);
        parcat(str, chi_angle[i][j], fmt, bstr);
        fprintf(fp, "%s\n", str);
    }
}

static void print_duplex_sugar_conformation(long num_bp, char **bp_seq,
                                            double **sugar_angle, char *bstr, char *fmt, FILE * fp)
{
    char str[BUF512];
    long i, j, k, ioffset, ds = 2;

    for (i = 1; i <= ds; i++) {
        (i == 1) ? fprintf(fp, "Strand I\n") : fprintf(fp, "Strand II\n");
        fprintf(fp, " base       v0      v1      v2      v3      v4"
                "      tm       P    Puckering\n");

        for (j = 1; j <= num_bp; j++) {
            sprintf(str, "%4ld %c ", j, bp_seq[i][j]);
            ioffset = (j - 1) * 7;
            for (k = 1; k <= 7; k++)
                parcat(str, sugar_angle[i][ioffset + k], fmt, bstr);
            add_sugar_pucker(sugar_angle[i][ioffset + 7], bstr, str);
            fprintf(fp, "%s\n", str);
        }

        if (i == 1)
            fprintf(fp, "\n");
    }
}

static void print_ss_sugar_conformation(long num_bp, char **bp_seq,
                                        double **sugar_angle, char *bstr, char *fmt, FILE * fp)
{
    char str[BUF512];
    long i = 1, j, k, ioffset;

    fprintf(fp, " base       v0      v1      v2      v3      v4"
            "      tm       P    Puckering\n");

    for (j = 1; j <= num_bp; j++) {
        sprintf(str, "%4ld %c ", j, bp_seq[i][j]);

        ioffset = (j - 1) * 7;
        for (k = 1; k <= 7; k++)
            parcat(str, sugar_angle[i][ioffset + k], fmt, bstr);
        add_sugar_pucker(sugar_angle[i][ioffset + 7], bstr, str);

        fprintf(fp, "%s\n", str);
    }
}

void backbone_torsion(long ds, long num_bp, long **pair_num, char **bp_seq, long **sugar,
                      long **chi, double **xyz, double **nt_bb_torsion, FILE * fp)
{
    char *bstr = "    --- ", *fmt = "%8.1f";
    double **chi_angle, **sugar_angle;
    long num_bpx7;

    num_bpx7 = num_bp * 7;

    chi_angle = dmatrix(1, ds, 1, num_bp);
    sugar_angle = dmatrix(1, ds, 1, num_bpx7);

    /* initialize with EMPTY_NUMBER */
    init_dmatrix(chi_angle, 1, ds, 1, num_bp, EMPTY_NUMBER);
    init_dmatrix(sugar_angle, 1, ds, 1, num_bpx7, EMPTY_NUMBER);
    get_chi_torsions(ds, num_bp, chi, xyz, chi_angle);
    get_sugar_torsions(ds, num_bp, sugar, xyz, sugar_angle);

    print_sep(fp, '*', 76);
    if (ds == 2)
        print_duplex_torsions(num_bp, pair_num, bp_seq, chi_angle, nt_bb_torsion, bstr, fmt, fp);
    else
        print_ss_torsions(num_bp, pair_num, bp_seq, chi_angle, nt_bb_torsion, bstr, fmt, fp);

    print_sep(fp, '*', 76);
    output_header_sugar_torsion(fp);
    if (ds == 2)
        print_duplex_sugar_conformation(num_bp, bp_seq, sugar_angle, bstr, fmt, fp);
    else
        print_ss_sugar_conformation(num_bp, bp_seq, sugar_angle, bstr, fmt, fp);

    free_dmatrix(chi_angle, 1, ds, 1, num_bp);
    free_dmatrix(sugar_angle, 1, ds, 1, num_bpx7);
}

/* Calculate same strand P-P, C1'-C1' distances */
void p_c1_dist(long ds, long num_bp, char **bp_seq, long **phos, long **chi,
               double **xyz, long *bphlx, FILE * fp)
{
    char *bstr = "       ---", *fmt = "%10.2f";
    char str[BUF512];
    double **c1_dist, **p_dist;
    long i, ia, ib, j, nbpm1;

    nbpm1 = num_bp - 1;

    p_dist = dmatrix(1, ds, 1, nbpm1);
    c1_dist = dmatrix(1, ds, 1, nbpm1);

    init_dmatrix(p_dist, 1, ds, 1, nbpm1, EMPTY_NUMBER);
    init_dmatrix(c1_dist, 1, ds, 1, nbpm1, EMPTY_NUMBER);

    for (i = 1; i <= ds; i++) {
        for (j = 1; j <= nbpm1; j++) {
            if (bphlx[j])  /* helix break */
                continue;
            ia = phos[i][j];
            ib = phos[i][j + 1];
            if (ia && ib)
                p_dist[i][j] = p1p2_dist(xyz[ia], xyz[ib]);
            ia = chi[i][(j - 1) * 4 + 2];
            ib = chi[i][j * 4 + 2];
            if (ia && ib)
                c1_dist[i][j] = p1p2_dist(xyz[ia], xyz[ib]);
        }
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "Same strand P--P and C1'--C1' virtual bond distances\n\n");

    fprintf(fp, "                 Strand I");
    if (ds == 2)
        fprintf(fp, "                          Strand II");
    fprintf(fp, "\n");

    fprintf(fp, "    step      P--P     C1'--C1'");
    if (ds == 2)
        fprintf(fp, "       step      P--P     C1'--C1'");
    fprintf(fp, "\n");

    for (i = 1; i <= nbpm1; i++) {
        j = 1;  /* strand I */
        sprintf(str, "%4ld %c/%c", i, bp_seq[j][i], bp_seq[j][i + 1]);
        parcat(str, p_dist[j][i], fmt, bstr);
        parcat(str, c1_dist[j][i], fmt, bstr);
        fprintf(fp, "%s", str);

        if (ds == 2) {
            sprintf(str, "      %4ld %c/%c", i, bp_seq[ds][i], bp_seq[ds][i + 1]);
            parcat(str, p_dist[ds][i], fmt, bstr);
            parcat(str, c1_dist[ds][i], fmt, bstr);
            fprintf(fp, "%s", str);
        }
        fprintf(fp, "\n");
    }

    free_dmatrix(p_dist, 1, ds, 1, nbpm1);
    free_dmatrix(c1_dist, 1, ds, 1, nbpm1);
}

/* Get lambda angle and C1'-C1', C6-C8, N1-N9 distances */
void lambda_d3(long num_bp, char **bp_seq, long **chi, long **c6_c8, double **xyz, FILE * fp)
{
    char str[BUF512], *bstr = "       ---", *fmt = "%10.1f";
    long i, j, ioffset, **c1_c1, **n1_n9;
    double vcc1[4], vcc2[4], vcn1[4], vcn2[4], **lambda_dist;

    c1_c1 = lmatrix(1, 2, 1, num_bp);
    n1_n9 = lmatrix(1, 2, 1, num_bp);
    lambda_dist = dmatrix(1, num_bp, 1, 5);

    init_dmatrix(lambda_dist, 1, num_bp, 1, 5, EMPTY_NUMBER);

    for (i = 1; i <= num_bp; i++) {
        ioffset = (i - 1) * 4;
        c1_c1[1][i] = chi[1][ioffset + 2];
        c1_c1[2][i] = chi[2][ioffset + 2];
        n1_n9[1][i] = chi[1][ioffset + 3];
        n1_n9[2][i] = chi[2][ioffset + 3];

        if (c1_c1[1][i] && c1_c1[2][i]) {
            ddxyz(xyz[c1_c1[2][i]], xyz[c1_c1[1][i]], vcc1);
            for (j = 1; j <= 3; j++)
                vcc2[j] = -vcc1[j];
            lambda_dist[i][3] = veclen(vcc1);  /* C1'-C1' distance */

            if (n1_n9[1][i] && n1_n9[2][i]) {
                ddxyz(xyz[c1_c1[1][i]], xyz[n1_n9[1][i]], vcn1);
                ddxyz(xyz[c1_c1[2][i]], xyz[n1_n9[2][i]], vcn2);
                lambda_dist[i][1] = magang(vcc2, vcn1);  /* lambda1 */
                lambda_dist[i][2] = magang(vcc1, vcn2);  /* lambda2 */

            } else if (n1_n9[1][i] && !n1_n9[2][i]) {
                ddxyz(xyz[c1_c1[1][i]], xyz[n1_n9[1][i]], vcn1);
                lambda_dist[i][1] = magang(vcc2, vcn1);  /* lambda1 */

            } else if (!n1_n9[1][i] && n1_n9[2][i]) {
                ddxyz(xyz[c1_c1[2][i]], xyz[n1_n9[2][i]], vcn2);
                lambda_dist[i][2] = magang(vcc1, vcn2);  /* lambda2 */
            }
        }

        if (n1_n9[1][i] && n1_n9[2][i])  /* N1-N9 distance */
            lambda_dist[i][4] = p1p2_dist(xyz[n1_n9[1][i]], xyz[n1_n9[2][i]]);
        if (c6_c8[1][i] && c6_c8[2][i])  /* C6-C8 distance */
            lambda_dist[i][5] = p1p2_dist(xyz[c6_c8[1][i]], xyz[c6_c8[2][i]]);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "lambda: virtual angle between C1'-YN1 or C1'-RN9"
            " glycosidic bonds and the\n"
            "        base-pair C1'-C1' line\n\n"
            "C1'-C1': distance between C1' atoms for each base-pair\n"
            "RN9-YN1: distance between RN9-YN1 atoms for each base-pair\n"
            "RC8-YC6: distance between RC8-YC6 atoms for each base-pair\n");

    fprintf(fp, "\n    bp     lambda(I) lambda(II)  C1'-C1'   RN9-YN1   RC8-YC6\n");
    for (i = 1; i <= num_bp; i++) {
        sprintf(str, "%4ld %c%c%c", i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
        for (j = 1; j <= 5; j++)
            parcat(str, lambda_dist[i][j], fmt, bstr);
        fprintf(fp, "%s\n", str);
    }

    /* write C1', RN9/YN1, RC8/YC6 xyz coordinates to "auxiliary.par" */
    print_axyz(num_bp, bp_seq, c1_c1, "C1'", xyz);
    print_axyz(num_bp, bp_seq, n1_n9, "RN9/YN1", xyz);
    print_axyz(num_bp, bp_seq, c6_c8, "RC8/YC6", xyz);

    free_lmatrix(c1_c1, 1, 2, 1, num_bp);
    free_lmatrix(n1_n9, 1, 2, 1, num_bp);
    free_dmatrix(lambda_dist, 1, num_bp, 1, 5);
}

/* Print xyz coordinates of P, C1', RN9/YN1 and RC8/YC6 */
void print_axyz(long num_bp, char **bp_seq, long **aidx, char *aname, double **xyz)
{
    char *bstr = "    ---- ", *fmt = "%9.3f";
    long i, j;
    FILE *fc;

    fc = open_file(AUX_FILE, "a");
    print_sep(fc, '*', 76);
    fprintf(fc, "xyz coordinates of %s atoms\n\n", aname);

    fprintf(fc, "    bp        xI       yI       zI       xII      yII      zII\n");
    for (i = 1; i <= num_bp; i++) {
        fprintf(fc, "%4ld %c%c%c ", i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
        if (aidx[1][i])
            for (j = 1; j <= 3; j++)
                fprintf(fc, fmt, xyz[aidx[1][i]][j]);
        else
            fprintf(fc, "%s%s%s", bstr, bstr, bstr);

        if (aidx[2][i])
            for (j = 1; j <= 3; j++)
                fprintf(fc, fmt, xyz[aidx[2][i]][j]);
        else
            fprintf(fc, "%s%s%s", bstr, bstr, bstr);

        fprintf(fc, "\n");
    }

    close_file(fc);
}

/* Calculate the correction angle for refined groove width */
static double gw_angle(long *idx, double *pvec12, double **xyz)
{
    double pvec1[4], pvec2[4], pvecm[4];

    ddxyz(xyz[idx[2]], xyz[idx[1]], pvec1);  /* strand I */
    ddxyz(xyz[idx[4]], xyz[idx[3]], pvec2);  /* strand II */
    vec_norm(pvec1);
    vec_norm(pvec2);
    sumxyz(pvec1, pvec2, pvecm);

    return magang(pvecm, pvec12);
}

/* Groove width parameters based on El Hassan and Calladine (1998) */
void groove_width(long parallel, long num_bp, char **bp_seq, long **phos,
                  double **xyz, long *bphlx, FILE * fp)
{
    char str[BUF512];
    char *bstr = "       ---", *fmt = "%10.1f";  /* groove width */
    char *bstr2 = "   ----", *fmt2 = "%7.1f";  /* P-P distance matrix */

    double anga, angb, dpa, dpb;
    double vecpa[4], vecpb[4];
    double **gwidth, **pdist;

    long iminor[5] = { BUF512, 1, -2, 2, -1 };
    long imajor[3] = { BUF512, -2, 2 };
    long idx[5], idxm1[5], idxp1[5], pp[5], ppa[5], ppb[5];
    long i, j, k, nbpm1, nset;
    long first_dist_num, items_per_line = 12, num_dist_sets;

    FILE *fc;

    nbpm1 = num_bp - 1;
    gwidth = dmatrix(1, nbpm1, 1, 4);
    init_dmatrix(gwidth, 1, nbpm1, 1, 4, EMPTY_NUMBER);

    /* method 1 is based on direct P-P distance
     *   minor: 0.5*[(P(i+1)-p(i-2))+(P(i+2)-p(i-1))]
     *   major: P(i-2)-p(i+2)
     * method 2 is method 1 plus a refinement
     *   minor: 0.5*[(P(i+1)-p(i-2))*sin(t1)+(P(i+2)-p(i-1))*sin(t2)]
     *   major: [P(i-2)-p(i+2)]*sin(t)
     */
    for (i = 1; i <= nbpm1; i++) {
        /* minor groove width */
        for (j = 1; j <= 4; j++) {
            idx[j] = i + iminor[j];
            if (!lval_in_range(idx[j], 1, nbpm1))
                break;
        }
        if (j > 4) {
            pp[1] = phos[1][idx[1] + 1];
            pp[3] = phos[1][idx[3] + 1];
            if (parallel) {
                pp[2] = phos[2][idx[2] + 1];
                pp[4] = phos[2][idx[4] + 1];
            } else {
                pp[2] = phos[2][idx[2]];
                pp[4] = phos[2][idx[4]];
            }

            if (pp[1] && pp[2] && pp[3] && pp[4]) {
                ddxyz(xyz[pp[2]], xyz[pp[1]], vecpa);
                ddxyz(xyz[pp[4]], xyz[pp[3]], vecpb);
                dpa = veclen(vecpa);
                dpb = veclen(vecpb);
                gwidth[i][1] = 0.5 * (dpa + dpb);  /* method 1 */

                /* method 2 (refined P-P distance) */
                for (k = 1; k <= 4; k++) {
                    idxm1[k] = idx[k] - 1;
                    idxp1[k] = idx[k] + 1;
                    if (idxm1[k] < 1 || idxp1[k] > nbpm1)
                        break;
                }
                if (k > 4) {
                    ppa[1] = phos[1][idxp1[1] + 1];
                    ppa[2] = phos[1][idxm1[1] + 1];
                    ppb[1] = phos[1][idxp1[3] + 1];
                    ppb[2] = phos[1][idxm1[3] + 1];
                    if (parallel) {
                        ppa[3] = phos[2][idxp1[2] + 1];
                        ppa[4] = phos[2][idxm1[2] + 1];
                        ppb[3] = phos[2][idxp1[4] + 1];
                        ppb[4] = phos[2][idxm1[4] + 1];
                    } else {
                        ppa[3] = phos[2][idxp1[2]];
                        ppa[4] = phos[2][idxm1[2]];
                        ppb[3] = phos[2][idxp1[4]];
                        ppb[4] = phos[2][idxm1[4]];
                    }

                    if (ppa[1] && ppa[2] && ppa[3] && ppa[4] &&
                        ppb[1] && ppb[2] && ppb[3] && ppb[4]) {
                        anga = deg2rad(gw_angle(ppa, vecpa, xyz));
                        angb = deg2rad(gw_angle(ppb, vecpb, xyz));
                        gwidth[i][2] = 0.5 * (dpa * sin(anga) + dpb * sin(angb));
                    }
                }
            }
        }

        /* major groove width */
        for (j = 1; j <= 2; j++) {
            idx[j] = i + imajor[j];
            if (!lval_in_range(idx[j], 1, nbpm1))
                break;
        }
        if (j > 2) {
            pp[1] = phos[1][idx[1] + 1];
            pp[2] = (parallel) ? phos[2][idx[2] + 1] : phos[2][idx[2]];
            if (pp[1] && pp[2]) {
                ddxyz(xyz[pp[2]], xyz[pp[1]], vecpa);
                dpa = veclen(vecpa);
                gwidth[i][3] = dpa;  /* method 1 */

                /* method 2 (refined P-P distance) */
                for (k = 1; k <= 2; k++) {
                    idxm1[k] = idx[k] - 1;
                    idxp1[k] = idx[k] + 1;
                    if (idxm1[k] < 1 || idxp1[k] > nbpm1)
                        break;
                }
                if (k > 2) {
                    ppa[1] = phos[1][idxp1[1] + 1];
                    ppa[2] = phos[1][idxm1[1] + 1];
                    if (parallel) {
                        ppa[3] = phos[2][idxp1[2] + 1];
                        ppa[4] = phos[2][idxm1[2] + 1];
                    } else {
                        ppa[3] = phos[2][idxp1[2]];
                        ppa[4] = phos[2][idxm1[2]];
                    }

                    if (ppa[1] && ppa[2] && ppa[3] && ppa[4]) {
                        anga = deg2rad(gw_angle(ppa, vecpa, xyz));
                        gwidth[i][4] = dpa * sin(anga);
                    }
                }
            }
        }
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "Minor and major groove widths: direct P-P distances "
            "and refined P-P distances\n   which take into account the "
            "directions of the sugar-phosphate backbones\n\n");
    fprintf(fp, "   (Subtract 5.8 Angstrom from the values to take account"
            " of the vdw radii\n    of the phosphate groups, and for"
            " comparison with FreeHelix and Curves.)\n\n");
    fprintf(fp, "Ref: M. A. El Hassan and C. R. Calladine (1998)."
            " ``Two Distinct Modes of\n     Protein-induced Bending"
            " in DNA.'' J. Mol. Biol., v282, pp331-343.\n\n");
    fprintf(fp, "                  Minor Groove        Major Groove\n"
            "                 P-P     Refined     P-P     Refined\n");

    /* taking into account helix breaks */
    for (i = 1; i <= nbpm1; i++) {
        if (bphlx[i]) {
            for (k = i - 3; k <= i + 3; k++) {
                if (!lval_in_range(k, 1, nbpm1))
                    continue;
                for (j = 1; j <= 4; j++)  /* for refined definition */
                    if (!((k == i - 3 || k == i + 3) && (j == 1 || j == 3)))
                        gwidth[k][j] = EMPTY_NUMBER;
            }
        }
    }

    for (i = 1; i <= nbpm1; i++) {
        sprintf(str, "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
                bp_seq[2][i + 1], bp_seq[2][i]);
        for (j = 1; j <= 4; j++)
            parcat(str, gwidth[i][j], fmt, bstr);
        fprintf(fp, "%s\n", str);
    }

    /* P-P distance matrix */
    pdist = dmatrix(1, num_bp, 1, num_bp);
    init_dmatrix(pdist, 1, num_bp, 1, num_bp, EMPTY_NUMBER);

    for (i = 1; i <= num_bp; i++)
        for (j = 1; j <= num_bp; j++)
            if (phos[1][i] && phos[2][j])
                pdist[i][j] = p1p2_dist(xyz[phos[1][i]], xyz[phos[2][j]]);
    fc = open_file(AUX_FILE, "a");

    print_sep(fc, '*', 91);
    fprintf(fc, "Phosphorus-phosphorus distance in Angstroms\n\n");
    num_dist_sets = (long) ceil(num_bp / (double) items_per_line);
    for (nset = 1; nset <= num_dist_sets; nset++) {
        if (nset == num_dist_sets) {
            k = num_bp % items_per_line;
            if (!k)
                k = items_per_line;
        } else
            k = items_per_line;

        first_dist_num = (nset - 1) * items_per_line;
        fprintf(fc, "      ");
        for (i = first_dist_num + 1; i <= first_dist_num + k; i++)
            fprintf(fc, "%7ld", i);
        fprintf(fc, "\n");

        fprintf(fc, "         ");
        for (i = first_dist_num + 1; i <= first_dist_num + k; i++)
            fprintf(fc, "   %c   ", bp_seq[2][i]);
        fprintf(fc, "\n");

        for (i = 1; i <= num_bp; i++) {
            sprintf(str, "%4ld %c ", i, bp_seq[1][i]);
            for (j = first_dist_num + 1; j <= first_dist_num + k; j++)
                parcat(str, pdist[i][j], fmt2, bstr2);
            fprintf(fc, "%s\n", str);
        }
        if (nset != num_dist_sets)
            fprintf(fc, "\n");
    }

    close_file(fc);

    free_dmatrix(gwidth, 1, nbpm1, 1, 4);
    free_dmatrix(pdist, 1, num_bp, 1, num_bp);
}

void check_wc_wobble_pair(long *bpid, char *bp, double shear, double stretch, double opening)
{
    static char *WC[9] = { WC_LIST };

    if (fabs(stretch) > 2.0 || fabs(opening) > 60)
        return;

    if (dval_in_range(fabs(shear), 1.8, 2.8))
        *bpid = 1;  /* with WC geometry */
    if (fabs(shear) <= 1.8 && num_strmatch(bp, WC, 1, 8))
        *bpid = 2;  /* WC */
}

/* Check if a base-pair is Watson-Crick:
   2: WC (1-below, plus |shear| <= 2.0, and base-pair sequence constraints)
   1: with correct WC geometry (i.e., x, y, z-axes in parallel directions)
   0: other cases, definitely non-WC (default) */
static void check_Watson_Crick(long num_bp, char **bp_seq, double **orien, double **org,
                               long *WC_info)
{
    char bpi[3];
    double pars[7], o1[4], o2[4], morg[4], **r1, **r2, **mst;
    long i, j, k;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);

    /* y- and z-axes of strand II base have been reversed */
    for (i = 1; i <= num_bp; i++) {
        j = (i - 1) * 9;
        sprintf(bpi, "%c%c", toupper((int) bp_seq[1][i]), toupper((int) bp_seq[2][i]));
        k = (bp_seq[0][i] == '-') &&  /* anti-parallel: y- & z-axes already REVERSED */
            dot(&orien[1][j], &orien[2][j]) > 0.0 &&  /* x-axis */
            dot(&orien[1][j + 3], &orien[2][j + 3]) > 0.0;  /* y-axis */
        if (k) {
            refs_right_left(i, orien, org, r1, o1, r2, o2);
            bpstep_par(r1, o1, r2, o2, pars, mst, morg);
            check_wc_wobble_pair(&WC_info[i], bpi, pars[1], pars[2], pars[6]);
        }
    }

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
}

/* see also set_nmarkers() in 'find_pair.c' */
void set_chain_nmarkers019_to_symbols(long num, long *nmarkers, char *cmarkers)
{
    char helix_begin = Gvars.CHAIN_MARKERS[0];
    char helix_middle = Gvars.CHAIN_MARKERS[1];
    char helix_end = Gvars.CHAIN_MARKERS[2];
    char isolated_bp = Gvars.CHAIN_MARKERS[3];
    long i, j, k, n, lastc;
    long *tcp, *temp, **seidx;

    lastc = strchr("0nofNOF", Gvars.CHAIN_MARKERS[4]) ? FALSE : TRUE;

    tcp = lvector(1, num);
    j = -1;  /* index for isolated bps */
    k = 1;  /* index for for helices */
    for (i = 1; i <= num; i++) {
        if (nmarkers[i] == 0)
            tcp[i] = k;
        else if (nmarkers[i] == 9) {
            tcp[i] = k;
            k++;
        } else {
            tcp[i] = j;
            j--;
        }
    }

    temp = lvector(1, num);
    for (i = 1; i < num; i++)
        temp[i] = (tcp[i + 1] != tcp[i]) ? 1 : 0;
    temp[num] = 1;

    n = 0;  /* get number of fragments */
    for (i = 1; i <= num; i++)
        if (temp[i])
            ++n;

    seidx = lmatrix(1, n, 0, 2);  /* allocate spaces */
    n = 0;
    for (i = 1; i <= num; i++)
        if (temp[i])
            seidx[++n][2] = i;
    for (i = 2; i <= n; i++)
        seidx[i][1] = seidx[i - 1][2] + 1;
    seidx[1][1] = 1;

    for (i = 1; i <= n; i++)
        seidx[i][0] = seidx[i][2] - seidx[i][1] + 1;

    for (i = 1; i <= n; i++) {
        k = seidx[i][0];
        if (k == 1) {  /* isolated bp */
            j = seidx[i][1];
            cmarkers[j] = isolated_bp;
        } else {
            for (j = seidx[i][1]; j <= seidx[i][2]; j++) {
                if (j == seidx[i][1])
                    cmarkers[j] = helix_begin;
                else if (j == seidx[i][2])
                    cmarkers[j] = helix_end;
                else
                    cmarkers[j] = helix_middle;
            }
        }
    }

    if (!lastc && cmarkers[seidx[n][2]] != isolated_bp)
        cmarkers[seidx[n][2]] = helix_middle;

    free_lvector(tcp, 1, num);
    free_lvector(temp, 1, num);
    free_lmatrix(seidx, 1, n, 0, 2);
}

void get_bp_3char_symbols(long bp_type, char zdir, char *bp_sym)
{
    sprintf(bp_sym, "%c%c%c", (bp_type == 2) ? '-' : '*', (bp_type > 0) ? '-' : '*', zdir);
}

/* Get the local reference frame for each base. Only the ring atoms are
   included in least-squares fitting */
void ref_frames(long ds, long num_bp, long **pair_num, char **bp_seq, long **seidx,
                long *RY, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                char **Miscs, double **xyz, FILE * fp, double **orien, double **org,
                long *WC_info, long *str_type, long irna, long **o3p_brk)
{
    static char *RingAtom[] = { RA_LIST };
    static char *rRingAtom[] =  /* Babcock also uses C1* atom! */
    { " C1'", RA_LIST };
    char bp_sym[BUF32], BDIR[BUF512], idmsg[BUF512], sidmsg[BUF512], spdb[BUF512];
    char *sChainID, **sAtomName, **sResName, **sMiscs, *cmarkers, *ss_markers;

    double orgi[4], vz[4];
    double **eRing_xyz, **fitted_xyz, **rms_fit, **sRing_xyz, **sxyz, **R;

    long i, ib, ie, ik, j, k, m, rnum, RingAtom_num, RA_NUM;
    long exp_katom, ioffset3, ioffset9, nmatch, snum, std_katom;
    long *sResSeq;

    irna ? get_BDIR(BDIR, "rAtomic_A.pdb") : get_BDIR(BDIR, "Atomic_A.pdb");

    RA_NUM = sizeof RingAtom / sizeof RingAtom[0];
    if (irna)
        RA_NUM++;  /* One more C1' atom */
    rms_fit = dmatrix(1, ds, 1, num_bp);
    eRing_xyz = dmatrix(1, RA_NUM, 1, 3);
    sRing_xyz = dmatrix(1, RA_NUM, 1, 3);
    sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
    sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
    sChainID = cvector(1, NUM_RESIDUE_ATOMS);
    sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
    sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
    sMiscs = cmatrix(1, NUM_RESIDUE_ATOMS, 0, NMISC);
    fitted_xyz = dmatrix(1, RA_NUM, 1, 3);
    R = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= ds; i++) {
        for (j = 1; j <= num_bp; j++) {
            rnum = pair_num[i][j];
            ib = seidx[rnum][1];
            ie = seidx[rnum][2];
            get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);
            RingAtom_num = (RY[rnum] == 1) ? RA_NUM : RA_NUM - 3;
            set_std_base_pdb(BDIR, irna, bp_seq[i][j], spdb);
            snum = read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz,
                            sMiscs, 1, "*");
            sprintf(sidmsg, "in standard base: %s", spdb);

            nmatch = 0;
            for (k = 0; k < RingAtom_num; k++) {
                if (irna) {
                    exp_katom = find_1st_atom(rRingAtom[k], AtomName, ib, ie, idmsg);
                    std_katom = find_1st_atom(rRingAtom[k], sAtomName, 1, snum, sidmsg);
                } else {
                    exp_katom = find_1st_atom(RingAtom[k], AtomName, ib, ie, idmsg);
                    std_katom = find_1st_atom(RingAtom[k], sAtomName, 1, snum, sidmsg);
                }
                if (exp_katom && std_katom) {
                    ++nmatch;
                    cpxyz(xyz[exp_katom], eRing_xyz[nmatch]);
                    cpxyz(sxyz[std_katom], sRing_xyz[nmatch]);
                }
            }

            rms_fit[i][j] = ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, orgi);

            ioffset9 = (j - 1) * 9;
            if (i == 2) {
                for (k = 1; k <= 3; k++)  /* column-wise */
                    vz[k] = R[k][3];
                if (dot(&orien[1][ioffset9 + 6], vz) < 0.0) {
                    bp_seq[0][j] = '-';  /* anti-parallel */
                    reverse_y_z_columns(R);
                } else
                    bp_seq[0][j] = '+';  /* parallel */
            }
            cpxyz(orgi, org[i] + (j - 1) * 3);
            mst2orien(orien[i], ioffset9, R);
        }
    }

    if (ds == 2) {
        check_Watson_Crick(num_bp, bp_seq, orien, org, WC_info);

        ib = 0;
        ik = 0;
        if (num_bp == 1)
            init_dvector(orgi, 1, 3, 0.0);
        for (i = 1; i <= ds; i++) {
            for (j = 1; j <= num_bp; j++) {
                ioffset3 = (j - 1) * 3;
                if (j < num_bp)
                    ddxyz(org[i] + ioffset3, org[i] + ioffset3 + 3, orgi);
                cpxyz(orien[i] + (j - 1) * 9 + 6, vz);
                if (!pair_num[3][j]) {
                    ik++;  /* non-breaks */
                    if (dot(orgi, vz) < 0.0 && WC_info[j])
                        ++ib;  /* z-axis reversed */
                }
            }
        }

        /* most likely left-handed Z-DNA */
        if (ib && ib == ik) {
            *str_type = 1;  /* with Z-axis reversed */
            for (i = 1; i <= ds; i++)
                for (j = 1; j <= num_bp; j++) {
                    ioffset9 = (j - 1) * 9;
                    negate_xyz(orien[i] + ioffset9);  /* reverse x-axis */
                    negate_xyz(orien[i] + ioffset9 + 6);  /* reverse z-axis */
                }
        }

        if (ib && ib != ik)
            *str_type = 2;  /* unusual cases */
        if (ik != ds * num_bp)
            *str_type = *str_type + 10;  /* more than one helices */
    }
    /* end of ds == 2 */

    /* write the least-squares fitting rms value */
    print_sep(fp, '*', 76);
    fprintf(fp, "RMSD of the bases");
    if (ds == 2)
        fprintf(fp, " (----- for WC bp, + for isolated bp, x for helix change)");

    fprintf(fp, "\n\n");
    fprintf(fp, "            Strand I");
    if (ds == 2)
        fprintf(fp, "                    Strand II          Helix");
    fprintf(fp, "\n");

    cmarkers = cvector(1, num_bp);
    ss_markers = cvector(1, num_bp);
    if (ds == 1)
        set_chain_nmarkers019_to_symbols(num_bp, o3p_brk[1], ss_markers);
    else
        set_chain_nmarkers019_to_symbols(num_bp, pair_num[ds + 1], cmarkers);

    for (i = 1; i <= num_bp; i++) {
        rnum = pair_num[1][i];
        ib = seidx[rnum][1];
        base_str(ChainID[ib], ResSeq[ib], Miscs[ib], ResName[ib], bp_seq[1][i], 1, idmsg);
        if (ds == 1)
            fprintf(fp, "%4ld   (%5.3f) %s     %c\n", i, rms_fit[1][i], idmsg, ss_markers[i]);
        else {
            rnum = pair_num[2][i];
            ib = seidx[rnum][1];
            base_str(ChainID[ib], ResSeq[ib], Miscs[ib], ResName[ib], bp_seq[2][i], 2, sidmsg);
            get_bp_3char_symbols(WC_info[i], bp_seq[0][i], bp_sym);
            fprintf(fp, "%4ld   (%5.3f) %s-%s-%s (%5.3f)     %c", i, rms_fit[1][i], idmsg,
                    bp_sym, sidmsg, rms_fit[2][i], cmarkers[i]);
            if (Gvars.VERBOSE)
                fprintf(fp, " %c-%c", o3p_brk[1][i] ? 'x' : '-', o3p_brk[2][i] ? 'x' : '-');
            fprintf(fp, "\n");
        }
    }

    if (ds == 2) {
        k = 0;
        m = 0;
        for (i = 1; i <= num_bp; i++) {
            if (WC_info[i] != 2)
                k++;
            if (!WC_info[i])
                m++;
        }
        if (k)
            fprintf(fp, "\nNote: This structure contains %ld[%ld] non-Watson-Crick"
                    " base-pair%s.\n", k, m, (k == 1) ? "" : "s");
    }

    free_dmatrix(rms_fit, 1, ds, 1, num_bp);
    free_dmatrix(eRing_xyz, 1, RA_NUM, 1, 3);
    free_dmatrix(sRing_xyz, 1, RA_NUM, 1, 3);
    free_pdb(NUM_RESIDUE_ATOMS, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, sMiscs);
    free_dmatrix(fitted_xyz, 1, RA_NUM, 1, 3);
    free_dmatrix(R, 1, 3, 1, 3);
    free_cvector(cmarkers, 1, num_bp);
    free_cvector(ss_markers, 1, num_bp);
}

/* Calculate step or base-pair parameters (using the CEHS scheme) */
void bpstep_par(double **rot1, double *org1, double **rot2, double *org2, double *pars,
                double **mst_orien, double *mst_org)
{
    double phi, rolltilt;
    double hinge[4], mstx[4], msty[4], mstz[4], t1[4], t2[4];
    double **para_bp1, **para_bp2, **temp;
    long i, j;

    for (i = 1; i <= 3; i++) {
        t1[i] = rot1[i][3];  /* z1 */
        t2[i] = rot2[i][3];  /* z2 */
    }

    cross(t1, t2, hinge);
    rolltilt = magang(t1, t2);

    /* Handle the special cases where t1 and t2 are perfectly
       anti-parallel, thus hinge can not be defined by t1 and t2. The
       vector would be (0 0 0), corresponding to tilt^2 + roll^2 = 180
       Date: July 20, 2005 */
    if (veclen(hinge) < XEPS && (fabs(rolltilt - 180) < XEPS || rolltilt < XEPS))
        for (i = 1; i <= 3; i++)
            hinge[i] = rot1[i][1] + rot2[i][1] + rot1[i][2] + rot2[i][2];

    para_bp1 = dmatrix(1, 3, 1, 3);
    para_bp2 = dmatrix(1, 3, 1, 3);
    temp = dmatrix(1, 3, 1, 3);

    arb_rotation(hinge, -0.5 * rolltilt, temp);
    multi_matrix(temp, 3, 3, rot2, 3, 3, para_bp2);
    arb_rotation(hinge, 0.5 * rolltilt, temp);
    multi_matrix(temp, 3, 3, rot1, 3, 3, para_bp1);

    for (i = 1; i <= 3; i++) {
        mstz[i] = para_bp2[i][3];  /* also para_bp1(:,3) */
        t1[i] = para_bp1[i][2];  /* y1 */
        t2[i] = para_bp2[i][2];  /* y2 */
    }

    /* twist is the angle between the two y- or x-axes */
    pars[6] = vec_ang(t1, t2, mstz);

    /* get y- and x-axes */
    get_vector(t1, mstz, 0.5 * pars[6], msty);
    cross(msty, mstz, mstx);

    avexyz(org1, org2, mst_org);
    ddxyz(org1, org2, t1);
    x_y_z_2_mtx(mstx, msty, mstz, mst_orien);

    /* get the xyz displacement parameters */
    for (i = 1; i <= 3; i++) {
        pars[i] = 0.0;
        for (j = 1; j <= 3; j++)
            pars[i] += t1[j] * mst_orien[j][i];
    }

    /* phi angle is defined by hinge and msty */
    phi = deg2rad(vec_ang(hinge, msty, mstz));

    /* get roll and tilt angles */
    pars[5] = rolltilt * cos(phi);
    pars[4] = rolltilt * sin(phi);

    free_dmatrix(para_bp1, 1, 3, 1, 3);
    free_dmatrix(para_bp2, 1, 3, 1, 3);
    free_dmatrix(temp, 1, 3, 1, 3);
}

/* Calculate local helical parameters & its middle frame */
void helical_par(double **rot1, double *org1, double **rot2, double *org2, double *pars,
                 double **mst_orien, double *mst_org)
{
    double AD_mag, phi, TipInc1, TipInc2, vlen;
    double axis_h[4], hinge1[4], hinge2[4], t1[4], t2[4];
    double AD_axis[4], org1_h[4], org2_h[4];
    double **rot1_h, **rot2_h, **temp;
    long i, j;

    for (i = 1; i <= 3; i++) {  /* column-wise */
        t1[i] = rot2[i][1] - rot1[i][1];  /* dx */
        t2[i] = rot2[i][2] - rot1[i][2];  /* dy */
    }
    cross(t1, t2, axis_h);
    vlen = veclen(axis_h);
    if (vlen < XEPS) {  /* for twist = 0.0 */
        axis_h[1] = 0.0;
        axis_h[2] = 0.0;
        axis_h[3] = 1.0;
    } else
        for (i = 1; i <= 3; i++)
            axis_h[i] /= vlen;

    temp = dmatrix(1, 3, 1, 3);

    rot1_h = dmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++)
        t1[i] = rot1[i][3];  /* z1 */
    TipInc1 = magang(axis_h, t1);
    cross(axis_h, t1, hinge1);
    arb_rotation(hinge1, -TipInc1, temp);
    multi_matrix(temp, 3, 3, rot1, 3, 3, rot1_h);

    rot2_h = dmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++)
        t2[i] = rot2[i][3];  /* z2 */
    TipInc2 = magang(axis_h, t2);
    cross(axis_h, t2, hinge2);
    arb_rotation(hinge2, -TipInc2, temp);
    multi_matrix(temp, 3, 3, rot2, 3, 3, rot2_h);

    for (i = 1; i <= 3; i++) {
        t1[i] = rot1_h[i][1] + rot2_h[i][1];  /* x1 + x2 */
        t2[i] = rot1_h[i][2] + rot2_h[i][2];  /* y1 + y2 */
    }
    vec_norm(t1);
    vec_norm(t2);
    x_y_z_2_mtx(t1, t2, axis_h, mst_orien);

    for (i = 1; i <= 3; i++) {
        t1[i] = rot1_h[i][2];  /* y1_h */
        t2[i] = rot2_h[i][2];  /* y2_h */
    }
    pars[6] = vec_ang(t1, t2, axis_h);

    ddxyz(org1, org2, t2);  /* org2-org1 */
    pars[3] = dot(t2, axis_h);

    phi = deg2rad(vec_ang(hinge1, t1, axis_h));
    pars[5] = TipInc1 * cos(phi);
    pars[4] = TipInc1 * sin(phi);

    for (i = 1; i <= 3; i++)
        t1[i] = t2[i] - pars[3] * axis_h[i];
    if (fabs(pars[6]) < HTWIST0)  /* twist = 0.0: cf <xhelfunc> */
        for (i = 1; i <= 3; i++)
            org1_h[i] = org1[i] + 0.5 * t1[i];
    else {
        get_vector(t1, axis_h, 90 - pars[6] / 2, AD_axis);
        AD_mag = 0.5 * veclen(t1) / sin(deg2rad(pars[6] / 2));
        for (i = 1; i <= 3; i++)
            org1_h[i] = org1[i] + AD_mag * AD_axis[i];
    }

    for (i = 1; i <= 3; i++)
        org2_h[i] = org1_h[i] + pars[3] * axis_h[i];
    avexyz(org1_h, org2_h, mst_org);
    ddxyz(org1_h, org1, t1);

    for (i = 1; i <= 2; i++) {
        pars[i] = 0.0;
        for (j = 1; j <= 3; j++)
            pars[i] += t1[j] * rot1_h[j][i];
    }

    free_dmatrix(rot1_h, 1, 3, 1, 3);
    free_dmatrix(temp, 1, 3, 1, 3);
    free_dmatrix(rot2_h, 1, 3, 1, 3);
}

/* Print base-pair, step and helical parameters */
void print_par(char **bp_seq, long num_bp, long ich, long ishel, double **param, FILE * fp)
{
    char *fmt = "%10.2f";
    double temp[7];
    long i, j;

    if (!lval_in_range(ich, 1, 4))
        fatal("wrong option [%ld] for printing parameters\n", ich);

    if (ich == 1) {  /* base-pair parameters */
        fprintf(fp, "     bp        Shear    Stretch   Stagger"
                "    Buckle  Propeller  Opening\n");
        for (i = 1; i <= num_bp; i++) {
            fprintf(fp, " %4ld %c%c%c ", i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
            for (j = 1; j <= 6; j++)
                fprintf(fp, fmt, param[i][j]);
            fprintf(fp, "\n");
        }

    } else {
        if (num_bp == 1)
            return;
        if (ishel)
            fprintf(fp, "    step       X-disp    Y-disp   h-Rise"
                    "     Incl.       Tip   h-Twist\n");
        else
            fprintf(fp, "    step       Shift     Slide      Rise"
                    "      Tilt      Roll     Twist\n");
        for (i = 1; i <= num_bp - 1; i++) {
            if (ich == 2)  /* for base-pair step */
                fprintf(fp, "%4ld %c%c/%c%c", i, bp_seq[1][i],
                        bp_seq[1][i + 1], bp_seq[2][i + 1], bp_seq[2][i]);
            else
                fprintf(fp, "%4ld  %c/%c ", i, bp_seq[ich - 2][i], bp_seq[ich - 2][i + 1]);
            for (j = 1; j <= 6; j++)
                fprintf(fp, fmt, param[i][j]);
            fprintf(fp, "\n");
        }
    }

    if (num_bp > 2) {
        j = (ich == 1) ? num_bp : num_bp - 1;
        fprintf(fp, "          ");
        print_sep(fp, '~', 60);
        fprintf(fp, "      ave.");
        ave_dmatrix(param, j, 6, temp);
        for (i = 1; i <= 6; i++)
            fprintf(fp, fmt, temp[i]);
        fprintf(fp, "\n");
        fprintf(fp, "      s.d.");
        std_dmatrix(param, j, 6, temp);
        for (i = 1; i <= 6; i++)
            fprintf(fp, fmt, temp[i]);
        fprintf(fp, "\n");
    }
}

/* Origin xyz coordinates followed by direction cosines of x-, y- & z-axes */
static void print_analyze_ref_frames(long ds, long num_bp, char **bp_seq, double *iorg,
                                     double *iorien, long **pair_num, char **nt_info)
{
    FILE *fp;
    long i, j, ia, ib, ik, ioffset3, ioffset9;

    fp = open_file(REF_FILE, "w");

    if (ds == 1)
        fprintf(fp, "%5ld bases\n", num_bp);
    else
        fprintf(fp, "%5ld base-pairs\n", num_bp);

    for (i = 1; i <= num_bp; i++) {
        ia = pair_num[1][i];
        if (ds == 1)
            fprintf(fp, "... %5ld %c   # %s\n", i, bp_seq[1][i], nt_info[ia]);
        else {
            ib = pair_num[2][i];
            fprintf(fp, "... %5ld %c%c%c   # %s - %s\n", i, bp_seq[1][i],
                    bp_seq[0][i], bp_seq[2][i], nt_info[ia], nt_info[ib]);
        }

        ioffset3 = (i - 1) * 3;
        fprintf(fp, "%10.4f %10.4f %10.4f  # origin\n", iorg[ioffset3 + 1],
                iorg[ioffset3 + 2], iorg[ioffset3 + 3]);
        ioffset9 = (i - 1) * 9;
        for (j = 1; j <= 3; j++) {
            ik = ioffset9 + (j - 1) * 3;
            fprintf(fp, "%10.4f %10.4f %10.4f  # %c-axis\n", iorien[ik + 1],
                    iorien[ik + 2], iorien[ik + 3], (j == 1) ? 'x' : (j == 2) ? 'y' : 'z');
        }
    }

    close_file(fp);
}

static void single_helix(long num_bp, char **bp_seq, double **step_par, double **heli_par,
                         double **orien, double **org, FILE * fp, long **pair_num, char **nt_info)
{
    char str[BUF512];
    double **parmtx;
    long nbpm1;
    FILE *fstep, *fheli, *fchek;

    nbpm1 = num_bp - 1;

    parmtx = dmatrix(1, nbpm1, 1, 6);

    print_sep(fp, '*', 76);
    fprintf(fp, "Local base step parameters\n");
    parvec2mtx(step_par[1], nbpm1, parmtx);
    print_par(bp_seq, num_bp, 3, 0, parmtx, fp);

    /* step parameters for rebuilding */
    fstep = open_file(BPSTEP_FILE, "w");
    sprintf(str, "%4ld # bases\n%4ld # ***local step parameters***\n", num_bp, 0L);
    strcat(str, "#      Shift     Slide     Rise      Tilt      Roll      Twist\n");
    print_ss_rebuild_pars(parmtx, num_bp, str, bp_seq, fstep);
    close_file(fstep);

    print_sep(fp, '*', 76);
    fprintf(fp, "Local base helical parameters\n");
    parvec2mtx(heli_par[1], nbpm1, parmtx);
    print_par(bp_seq, num_bp, 3, 1, parmtx, fp);

    /* helical parameters for rebuilding */
    fheli = open_file(HLXSTEP_FILE, "w");
    sprintf(str, "%4ld # bases\n%4ld # ***local helical parameters***\n", num_bp, 1L);
    strcat(str, "#      X-disp    Y-disp    h-Rise    Incl.     Tip     h-Twist\n");
    print_ss_rebuild_pars(parmtx, num_bp, str, bp_seq, fheli);
    close_file(fheli);

    /* for checking out */
    fchek = open_file(AUX_FILE, "w");
    fprintf(fchek, "Reference frame: Origins (Ox, Oy, Oz) followed by the"
            " direction cosines of the\n"
            "                 X- (Xx, Yx, Zx), Y- (Yx, Yx, Yx)," " and Z- (Zx, Zx, Zx) axes\n");
    print_sep(fchek, '*', 89);
    fprintf(fchek, "Local base reference frames\n");
    print_ref(bp_seq, num_bp, 2, org[1], orien[1], fchek);
    close_file(fchek);

    /* reference frame for reseting the structure */
    print_analyze_ref_frames(1, num_bp, bp_seq, org[1], orien[1], pair_num, nt_info);

    free_dmatrix(parmtx, 1, num_bp, 1, 6);
}

void output_ave_std(long num, double **parcln, int dnum, char *fmt, FILE * fp)
{
    double temp[7];
    int k = 10;
    long i;

    if (num <= 2)
        return;

    k += dnum;
    fprintf(fp, "%*s", k, " ");
    print_sep(fp, '~', 60);

    fprintf(fp, "%*s      ave.", dnum, " ");
    ave_dmatrix(parcln, num, 6, temp);
    for (i = 1; i <= 6; i++)
        fprintf(fp, fmt, temp[i]);
    fprintf(fp, "\n");

    fprintf(fp, "%*s      s.d.", dnum, " ");
    std_dmatrix(parcln, num, 6, temp);
    for (i = 1; i <= 6; i++)
        fprintf(fp, fmt, temp[i]);
    fprintf(fp, "\n");
}

/* Print local base-pair step and helical parameter with helix breaks deleted */
void prt_stepstr(char **step_str, long num_step, long *bphlx, long ishel, double **param,
                 FILE * fp)
{
    char *bstr = "      ----", *fmt = "%10.2f";
    double **parcln;
    long i, j, num = 0;

    if (!num_step)  /* empty */
        return;

    parcln = dmatrix(1, num_step, 1, 6);
    if (ishel)
        fprintf(fp, "    step       X-disp    Y-disp   h-Rise" "     Incl.       Tip   h-Twist\n");
    else
        fprintf(fp, "    step       Shift     Slide      Rise" "      Tilt      Roll     Twist\n");
    for (i = 1; i <= num_step; i++) {
        fprintf(fp, "%4ld %s", i, step_str[i]);
        if (bphlx[i])
            for (j = 1; j <= 6; j++)
                fprintf(fp, "%s", bstr);
        else {
            num++;
            for (j = 1; j <= 6; j++) {
                fprintf(fp, fmt, param[i][j]);
                parcln[num][j] = param[i][j];
            }
        }
        fprintf(fp, "\n");
    }

    output_ave_std(num, parcln, 0, fmt, fp);

    free_dmatrix(parcln, 1, num_step, 1, 6);
}

/* Print local base-pair step and helical parameter with helix breaks
 * deleted. Some heavily kinked steps (e.g., pdt040) could be divided */
void prt_step_par(char **bp_seq, long num_bp, long *bphlx, long ishel, double **param, FILE * fp)
{
    char *bstr = "      ----", *fmt = "%10.2f";
    double temp[7], **parcln;
    long i, j, num = 0, nbpm1;

    if (num_bp == 1)
        return;

    nbpm1 = num_bp - 1;
    parcln = dmatrix(1, nbpm1, 1, 6);
    if (ishel)
        fprintf(fp, "    step       X-disp    Y-disp   h-Rise" "     Incl.       Tip   h-Twist\n");
    else
        fprintf(fp, "    step       Shift     Slide      Rise" "      Tilt      Roll     Twist\n");
    for (i = 1; i <= nbpm1; i++) {
        fprintf(fp, "%4ld %c%c/%c%c", i, bp_seq[1][i],
                bp_seq[1][i + 1], bp_seq[2][i + 1], bp_seq[2][i]);
        if (bphlx[i])
            for (j = 1; j <= 6; j++)
                fprintf(fp, "%s", bstr);
        else {
            num++;
            for (j = 1; j <= 6; j++) {
                fprintf(fp, fmt, param[i][j]);
                parcln[num][j] = param[i][j];
            }
        }
        fprintf(fp, "\n");
    }

    if (num > 2) {
        fprintf(fp, "          ");
        print_sep(fp, '~', 60);
        fprintf(fp, "      ave.");
        ave_dmatrix(parcln, num, 6, temp);
        for (i = 1; i <= 6; i++)
            fprintf(fp, fmt, temp[i]);
        fprintf(fp, "\n");
        fprintf(fp, "      s.d.");
        std_dmatrix(parcln, num, 6, temp);
        for (i = 1; i <= 6; i++)
            fprintf(fp, fmt, temp[i]);
        fprintf(fp, "\n");
    }
    free_dmatrix(parcln, 1, nbpm1, 1, 6);
}

static void double_helix(long num_bp, char **bp_seq, double **step_par, double **heli_par,
                         double **orien, double **org, long *WC_info, FILE * fp,
                         double **twist_rise, double *mst_orien, double *mst_org,
                         double *mst_orienH, double *mst_orgH, long *bphlx, long istart,
                         long istep, long bz, long *str_type, long **pair_num, char **nt_info)
{
    char str[BUF512], *fmt = " %9.3f";
    char **step_str;
    double hfoi[4], mfoi[4], o1[4], o2[4];
    double *bp_org, *bp_orien, **bp_par, **bp_step_par, **bp_heli_par;
    double **mfi, **hfi, **r1, **r2;
    long num_step = 0, bz_junction = 0, z_step = 0;
    long i, ib, ie, ioffset3, ioffset9, j, nbpm1;
    FILE *fchek, *fheli, *fstep;

    nbpm1 = num_bp - 1;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mfi = dmatrix(1, 3, 1, 3);
    hfi = dmatrix(1, 3, 1, 3);

    bp_par = dmatrix(1, num_bp, 1, 6);
    bp_org = dvector(1, num_bp * 3);
    bp_orien = dvector(1, num_bp * 9);

    print_sep(fp, '*', 76);
    fprintf(fp, "Origin (Ox, Oy, Oz) and mean normal vector"
            " (Nx, Ny, Nz) of each base-pair in\n"
            "   the coordinate system of the given structure\n\n");
    fprintf(fp, "      bp        Ox        Oy        Oz        Nx        Ny        Nz\n");

    for (i = 1; i <= num_bp; i++) {
        refs_right_left(i, orien, org, r1, o1, r2, o2);
        bpstep_par(r1, o1, r2, o2, bp_par[i], mfi, mfoi);

        fprintf(fp, " %4ld %c%c%c ", i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
        for (j = 1; j <= 3; j++)
            fprintf(fp, fmt, mfoi[j]);  /* origin */
        for (j = 1; j <= 3; j++)
            fprintf(fp, fmt, mfi[j][3]);  /* base-pair normal */
        fprintf(fp, "\n");

        cpxyz(mfoi, bp_org + (i - 1) * 3);
        mst2orien(bp_orien, (i - 1) * 9, mfi);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "Local base-pair parameters\n");
    print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

    bp_step_par = dmatrix(1, nbpm1, 1, 6);
    bp_heli_par = dmatrix(1, nbpm1, 1, 6);
    step_str = cmatrix(1, num_bp, 0, 5);
    ie = istart;
    for (i = istart; i <= nbpm1; i++) {
        if (istep > 0) {  /* continuous 1 to 3; 2 to 4 etc */
            ib = i;
            ie = ib + istep;
        } else {  /* next segment: 1 to 3, 3 to 5 etc */
            ib = ie;
            ie = ib - istep;
        }
        if (ie > num_bp)
            break;
        num_step++;
        sprintf(step_str[num_step], "%c%c/%c%c", bp_seq[1][ib],
                bp_seq[1][ie], bp_seq[2][ie], bp_seq[2][ib]);
        refs_i_j(ib, ie, bp_orien, bp_org, r1, o1, r2, o2);
        if (WC_info[i] == 2 && WC_info[i + 1] == 2)  /* must be WC step */
            bz_check(r1, o1, r2, o2, bz, &bz_junction, &z_step);
        bpstep_par(r1, o1, r2, o2, bp_step_par[num_step], mfi, mfoi);
        helical_par(r1, o1, r2, o2, bp_heli_par[num_step], hfi, hfoi);

        /* this could be meaningless for non-dinucleotide steps */
        twist_rise[num_step][1] = bp_step_par[num_step][6];
        twist_rise[num_step][2] = bp_step_par[num_step][3];

        ioffset3 = (num_step - 1) * 3;
        cpxyz(mfoi, mst_org + ioffset3);
        cpxyz(hfoi, mst_orgH + ioffset3);

        ioffset9 = (num_step - 1) * 9;
        mst2orien(mst_orien, ioffset9, mfi);
        mst2orien(mst_orienH, ioffset9, hfi);
    }

    if (bz && bz_junction) {
        *str_type = -2;  /* A-|B-/Z-DNA junction */
        if (z_step > 1)
            *str_type = -12;
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "Local base-pair step parameters\n");
    prt_stepstr(step_str, num_step, bphlx, 0, bp_step_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "Local base-pair helical parameters\n");
    prt_stepstr(step_str, num_step, bphlx, 1, bp_heli_par, fp);

    if (istep != 1) {
        print_sep(fp, '*', 76);
        fprintf(fp, "Step size [%ld]:  ", istep);
        fprintf(fp, "dinucleotide step classification based on Zp and ZpH will\n");
        fprintf(fp, "    be meaningless. Files <bp_helical.par>, <bp_step.par>, "
                "<hstacking.pdb>,\n");
        fprintf(fp, "    and <stacking.pdb> will also be corrupted.\n");
    }

    /* base-pair & step parameters for rebuilding */
    fstep = open_file(BPSTEP_FILE, "w");
    sprintf(str, "%4ld # base-pairs\n"
            "%4ld # ***local base-pair & step parameters***\n", num_bp, 0L);
    strcat(str, "#        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening"
           "     Shift     Slide     Rise      Tilt      Roll      Twist\n");
    print_ds_rebuild_pars(bp_par, bp_step_par, num_bp, str, bp_seq, fstep);
    close_file(fstep);

    /* base-pair & helical parameters for rebuilding */
    fheli = open_file(HLXSTEP_FILE, "w");
    sprintf(str, "%4ld # base-pairs\n"
            "%4ld # ***local base-pair & helical parameters***\n", num_bp, 1L);
    strcat(str, "#        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening"
           "     X-disp    Y-disp    h-Rise    Incl.     Tip     h-Twist\n");
    print_ds_rebuild_pars(bp_par, bp_heli_par, num_bp, str, bp_seq, fheli);
    close_file(fheli);

    /* for checking out */
    fchek = open_file(AUX_FILE, "w");
    fprintf(fchek, "Reference frame: Origins (Ox, Oy, Oz) followed by the"
            " direction cosines of the\n"
            "                 X- (Xx, Yx, Zx), Y- (Yx, Yx, Yx)," " and Z- (Zx, Zx, Zx) axes\n");

    print_sep(fchek, '*', 89);
    fprintf(fchek, "Local base-pair reference frames\n");
    print_ref(bp_seq, num_bp, 1, bp_org, bp_orien, fchek);

    print_sep(fchek, '*', 89);
    fprintf(fchek, "Local middle reference frames\n");
    print_ref(bp_seq, num_bp - 1, 4, mst_org, mst_orien, fchek);

    print_sep(fchek, '*', 89);
    fprintf(fchek, "Local middle helical reference frames\n");
    print_ref(bp_seq, num_bp - 1, 4, mst_orgH, mst_orienH, fchek);

    print_sep(fchek, '*', 89);
    fprintf(fchek, "Local strand I base reference frames\n");
    print_ref(bp_seq, num_bp, 2, org[1], orien[1], fchek);

    print_sep(fchek, '*', 89);
    fprintf(fchek, "Local strand II base reference frames\n");
    print_ref(bp_seq, num_bp, 3, org[2], orien[2], fchek);

    print_sep(fchek, '*', 76);
    fprintf(fchek, "Local strand I base step parameters\n");
    parvec2mtx(step_par[1], nbpm1, bp_step_par);
    print_par(bp_seq, num_bp, 3, 0, bp_step_par, fchek);

    print_sep(fchek, '*', 76);
    fprintf(fchek, "Local strand I base helical parameters\n");
    parvec2mtx(heli_par[1], nbpm1, bp_heli_par);
    print_par(bp_seq, num_bp, 3, 1, bp_heli_par, fchek);

    print_sep(fchek, '*', 76);
    fprintf(fchek, "Local strand II base step parameters\n");
    parvec2mtx(step_par[2], nbpm1, bp_step_par);
    print_par(bp_seq, num_bp, 4, 0, bp_step_par, fchek);

    print_sep(fchek, '*', 76);
    fprintf(fchek, "Local strand II base helical parameters\n");
    parvec2mtx(heli_par[2], nbpm1, bp_heli_par);
    print_par(bp_seq, num_bp, 4, 1, bp_heli_par, fchek);

    close_file(fchek);

    /* reference frame for reseting the structure */
    print_analyze_ref_frames(2, num_bp, bp_seq, bp_org, bp_orien, pair_num, nt_info);

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mfi, 1, 3, 1, 3);
    free_dmatrix(hfi, 1, 3, 1, 3);
    free_dmatrix(bp_step_par, 1, nbpm1, 1, 6);
    free_dmatrix(bp_heli_par, 1, nbpm1, 1, 6);
    free_dmatrix(bp_par, 1, num_bp, 1, 6);
    free_dvector(bp_org, 1, num_bp * 3);
    free_dvector(bp_orien, 1, num_bp * 9);
    free_cmatrix(step_str, 1, num_bp, 0, 5);
}

void bz_check(double **r1, double *o1, double **r2, double *o2, long bz,
              long *bz_junction, long *z_step)
{
    double x1[4], y1_osx[4], z1[4], x2[4], y2[4], z2[4], dorg[4];
    long i;

    if (!bz)  /* without bz option */
        return;

    for (i = 1; i <= 3; i++) {
        x1[i] = r1[i][1];
        y1_osx[i] = r1[i][2];
        z1[i] = r1[i][3];
        x2[i] = r2[i][1];
        y2[i] = r2[i][2];
        z2[i] = r2[i][3];
        dorg[i] = o2[i] - o1[i];
    }

    if (dot(x1, x2) < 0.0 && dot(z1, z2) < 0.0 && dot(y1_osx, y2) > 0.0) {
        if (dot(dorg, z1) > 0.0) {  /* bp2 in Z-form */
            for (i = 1; i <= 3; i++) {
                r2[i][1] = -r2[i][1];
                r2[i][3] = -r2[i][3];
            }
        } else {  /* bp1 in Z-form */
            for (i = 1; i <= 3; i++) {
                r1[i][1] = -r1[i][1];
                r1[i][3] = -r1[i][3];
            }
        }
        (*bz_junction)++;
        return;
    }

    if (dot(x1, x2) > 0.0 && dot(z1, z2) > 0.0 && dot(y1_osx, y2) > 0.0 && dot(dorg, z1) < 0.0 && dot(dorg, z2) < 0.0) {  /* Z-step */
        for (i = 1; i <= 3; i++) {
            r1[i][1] = -r1[i][1];  /* x1 */
            r1[i][3] = -r1[i][3];  /* z1 */
            r2[i][1] = -r2[i][1];  /* x2 */
            r2[i][3] = -r2[i][3];  /* z2 */
        }
        (*z_step)++;
    }
}

static double chk_twist_for_larger_segment(double sum, long num, long nbpm1, long *idx)
{
    long i;
    double ave = 0.0;

    if (num) {
        ave = sum / num;
        for (i = 1; i < nbpm1; i++)
            if (idx[i] && idx[i + 1])
                break;
        if (i == nbpm1 && i != 1)
            ave = 0.0;
    }

    return ave;
}

/* Get mean base-pair step twist angle: excluding breaks and non_WC steps */
void get_mtwist(long nbpm1, long *bphlx, long *WC_info, double **twist_rise,
                double *twist_p, double *twist_n)
{
    long i, num_p = 0, num_n = 0, *idx_p, *idx_n;
    double twist, tp = 0.0, tn = 0.0;

    idx_p = lvector(1, nbpm1);
    idx_n = lvector(1, nbpm1);

    for (i = 1; i <= nbpm1; i++) {
        if (!bphlx[i] && WC_info[i] == 2 && WC_info[i + 1] == 2) {
            twist = twist_rise[i][1];
            if (twist >= 0) {
                idx_p[i] = 1;
                num_p++;
                tp += twist;
            } else {
                idx_n[i] = 1;
                num_n++;
                tn += twist;
            }
        }
    }

    *twist_p = chk_twist_for_larger_segment(tp, num_p, nbpm1, idx_p);
    *twist_n = chk_twist_for_larger_segment(tn, num_n, nbpm1, idx_n);

    free_lvector(idx_p, 1, DUMMY);
    free_lvector(idx_n, 1, DUMMY);
}

/* Calculate and print out 3DNA recommended local parameters */
void get_parameters(long ds, long num_bp, char **bp_seq, double **orien, double **org,
                    long *WC_info, FILE * fp, double **twist_rise, double *mst_orien,
                    double *mst_org, double *mst_orienH, double *mst_orgH, long *bphlx,
                    long istart, long istep, long bz, long *str_type, long **pair_num,
                    char **nt_info)
{
    double hfoi[4], mfoi[4], o1[4], o2[4];
    double **heli_par, **step_par;
    double **hfi, **mfi, **r1, **r2;
    long i, j, m, ioffset3, ioffset9, nbpm1;

    nbpm1 = num_bp - 1;

    /* step and helical parameters for each strand */
    step_par = dmatrix(1, ds, 1, nbpm1 * 6);
    heli_par = dmatrix(1, ds, 1, nbpm1 * 6);
    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mfi = dmatrix(1, 3, 1, 3);
    hfi = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= ds; i++) {
        for (j = 1; j <= nbpm1; j++) {
            refs_i_j(j, j + 1, orien[i], org[i], r1, o1, r2, o2);
            m = (j - 1) * 6;
            bpstep_par(r1, o1, r2, o2, step_par[i] + m, mfi, mfoi);
            helical_par(r1, o1, r2, o2, heli_par[i] + m, hfi, hfoi);

            if (ds == 1) {  /* populate mst_orien/mst_org, mst_orienH/mst_orgH */
                twist_rise[j][1] = step_par[i][m + 6];
                twist_rise[j][2] = step_par[i][m + 3];

                ioffset3 = (j - 1) * 3;
                cpxyz(mfoi, mst_org + ioffset3);
                cpxyz(hfoi, mst_orgH + ioffset3);

                ioffset9 = (j - 1) * 9;
                mst2orien(mst_orien, ioffset9, mfi);
                mst2orien(mst_orienH, ioffset9, hfi);
            }
        }
    }

    if (ds == 1)
        single_helix(num_bp, bp_seq, step_par, heli_par, orien, org, fp, pair_num, nt_info);
    else
        double_helix(num_bp, bp_seq, step_par, heli_par, orien, org, WC_info, fp,
                     twist_rise, mst_orien, mst_org, mst_orienH, mst_orgH, bphlx,
                     istart, istep, bz, str_type, pair_num, nt_info);

    free_dmatrix(step_par, 1, ds, 1, nbpm1 * 6);
    free_dmatrix(heli_par, 1, ds, 1, nbpm1 * 6);
    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mfi, 1, 3, 1, 3);
    free_dmatrix(hfi, 1, 3, 1, 3);
}

/* Change vector-wise parameters to a num-by-6 matrix */
void parvec2mtx(double *parvec, long num, double **parmtx)
{
    long i, ioffset, j, nc = 6;

    for (i = 1; i <= num; i++) {
        ioffset = (i - 1) * nc;
        for (j = 1; j <= nc; j++)
            parmtx[i][j] = parvec[ioffset + j];
    }
}

/* Print parameters for the rebuilding of a single helix */
void print_ss_rebuild_pars(double **pars, long num_bp, char *str, char **bp_seq, FILE * fp)
{
    char *fmt = " %9.3f";
    long i, j;

    fprintf(fp, "%s", str);

    /* 1st base: 6 zeros */
    fprintf(fp, "%c ", bp_seq[1][1]);
    for (i = 1; i <= 6; i++)
        fprintf(fp, fmt, 0.0);
    fprintf(fp, "\n");

    for (i = 2; i <= num_bp; i++) {
        fprintf(fp, "%c ", bp_seq[1][i]);
        for (j = 1; j <= 6; j++)
            fprintf(fp, fmt, pars[i - 1][j]);
        fprintf(fp, "\n");
    }
}

/* Print parameters for the rebuilding of a duplex */
void print_ds_rebuild_pars(double **bp_par, double **step_par, long num_bp, char *str,
                           char **bp_seq, FILE * fp)
{
    char *fmt = " %9.3f";
    long i, j;

    fprintf(fp, "%s", str);

    /* 1st base-pair: 6 zeros for step parameters */
    fprintf(fp, "%c%c%c ", bp_seq[1][1], bp_seq[0][1], bp_seq[2][1]);
    for (i = 1; i <= 6; i++)
        fprintf(fp, fmt, bp_par[1][i]);
    for (i = 1; i <= 6; i++)
        fprintf(fp, fmt, 0.0);
    fprintf(fp, "\n");

    for (i = 2; i <= num_bp; i++) {
        fprintf(fp, "%c%c%c ", bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
        for (j = 1; j <= 6; j++)
            fprintf(fp, fmt, bp_par[i][j]);
        for (j = 1; j <= 6; j++)
            fprintf(fp, fmt, step_par[i - 1][j]);
        fprintf(fp, "\n");
    }
}

/* Print local base and base-pair reference frames */
void print_ref(char **bp_seq, long num_item, long ich, double *org, double *orien, FILE * fp)
{
    long i, ioffset3, ioffset9, j;

    if (!lval_in_range(ich, 1, 4))
        fatal("wrong option [%ld] for printing reference frames\n", ich);

    fprintf(fp, "                Ox      Oy      Oz     Xx    Xy    Xz"
            "   Yx    Yy    Yz    Zx    Zy    Zz\n");

    for (i = 1; i <= num_item; i++) {
        if (ich == 1)  /* base-pair */
            fprintf(fp, "%4ld %c%c%c   ", i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
        else if (ich == 4)  /* step */
            fprintf(fp, "%4ld %c%c/%c%c ", i, bp_seq[1][i],
                    bp_seq[1][i + 1], bp_seq[2][i + 1], bp_seq[2][i]);
        else  /* base I or II */
            fprintf(fp, "%4ld %c     ", i, bp_seq[ich - 1][i]);

        ioffset3 = (i - 1) * 3;
        for (j = 1; j <= 3; j++)
            fprintf(fp, "%8.2f", org[ioffset3 + j]);

        ioffset9 = (i - 1) * 9;
        for (j = 1; j <= 9; j++)
            fprintf(fp, "%6.2f", orien[ioffset9 + j]);

        fprintf(fp, "\n");
    }
}

/* Write multiple dinucleotide structures w.r.t. middle frames */
void write_mst(long ds, long num_bp, long **pair_num, char **bp_seq, double *mst_orien,
               double *mst_org, long **seidx, char **AtomName, char **ResName,
               char *ChainID, long *ResSeq, double **xyz, char **Miscs,
               long **htm_water, double **twist_rise, char *strfile)
{
    double rise, **mst, **xyz_residue;
    long i, inum, ioffset3, j, jr, k, m;
    long tnum_res1, tnum_res2, inum_base;
    long ivec[3], ivec0[BUF512], ivect[BUF512];
    FILE *fp;

    if (num_bp == 1)
        return;

    fp = open_file(strfile, "w");

    xyz_residue = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
    mst = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= num_bp - 1; i++) {
        rise = twist_rise[i][2];
        inum = 0;
        ioffset3 = (i - 1) * 3;
        orien2mst(mst_orien, (i - 1) * 9, mst);

        fprintf(fp, "%6s    %4ld\n", "MODEL ", i);
        if (!dval_in_range(rise, 2.0, 9.0))  /* rise out of normal range */
            fprintf(fp, "REMARK    NB [%.2f] -- the following step is unlikely in"
                    " stacking geometry!\n", rise);

        if (ds == 1)
            fprintf(fp, "REMARK    Section #%4.4ld %c/%c\n", i, bp_seq[1][i], bp_seq[1][i + 1]);
        else
            fprintf(fp, "REMARK    Section #%4.4ld %c%c/%c%c\n", i,
                    bp_seq[1][i], bp_seq[1][i + 1], bp_seq[2][i + 1], bp_seq[2][i]);

        fprintf(fp, "REMARK    %s\n", Gvars.X3DNA_VER);

        ivec[1] = pair_num[1][i];  /* lower level */
        if (ds == 1)
            inum_base = 1;
        else {
            inum_base = 2;
            ivec[2] = pair_num[2][i];
        }
        tnum_res1 = attached_residues(inum_base, ivec, ivec0, seidx, xyz, htm_water,
                                      &Gvars.misc_pars);
        fprintf(fp, "REMARK    LOWER: %4ld", tnum_res1);
        for (j = 1; j <= tnum_res1; j++) {
            ivect[j] = ivec0[j];
            fprintf(fp, "%4ld", j);
        }
        fprintf(fp, "\n");

        ivec[1] = pair_num[1][i + 1];  /* upper level */
        if (ds == 1)
            inum_base = 1;
        else {
            inum_base = 2;
            ivec[2] = pair_num[2][i + 1];
        }
        tnum_res2 = attached_residues(inum_base, ivec, ivec0, seidx, xyz, htm_water,
                                      &Gvars.misc_pars);
        fprintf(fp, "REMARK    UPPER: %4ld", tnum_res2);
        for (j = 1; j <= tnum_res2; j++) {
            k = j + tnum_res1;
            ivect[k] = ivec0[j];
            fprintf(fp, "%4ld", k);
        }
        fprintf(fp, "\n");

        for (j = 1; j <= tnum_res1 + tnum_res2; j++) {
            jr = ivect[j];
            for (k = seidx[jr][1]; k <= seidx[jr][2]; k++) {
                m = k - seidx[jr][1] + 1;
                cpxyz(xyz[k], xyz_residue[m]);
            }
            change_xyz(0, &mst_org[ioffset3], mst, seidx[jr][2] - seidx[jr][1] + 1, xyz_residue);
            pdb_record(seidx[jr][1], seidx[jr][2], &inum, 1, AtomName, ResName,
                       ChainID, ResSeq, xyz_residue, Miscs, fp);
        }
        fprintf(fp, "ENDMDL\n");
    }

    free_dmatrix(xyz_residue, 1, NUM_RESIDUE_ATOMS, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);

    close_file(fp);
}

void print_xyzP(long parallel, long nbpm1, char **bp_seq, long **phos, double *mst_orien,
                double *mst_org, double **xyz, FILE * fp, char *title_str, double **aveP,
                long p_offset)
{
    char *bstr = "    --- ", *fmt = "%8.2f%8.2f%8.2f";
    char str[BUF512], temp_str[BUF512];
    double P_mst1[4], P_mst2[4], temp[4];
    long i, ioffset3, ioffset9, ip1, ip2, j;

    print_sep(fp, '*', 76);
    fprintf(fp, "%s", title_str);

    for (i = 1; i <= nbpm1; i++) {
        sprintf(str, "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
                bp_seq[2][i + 1], bp_seq[2][i]);

        ioffset3 = (i - 1) * 3;
        ioffset9 = (i - 1) * 9;

        ip1 = phos[1 + p_offset][i + 1];  /* strand I */
        if (ip1) {
            ddxyz(mst_org + ioffset3, xyz[ip1], temp);
            for (j = 1; j <= 3; j++)
                P_mst1[j] = dot(temp, &mst_orien[ioffset9 + (j - 1) * 3]);
            sprintf(temp_str, fmt, P_mst1[1], P_mst1[2], P_mst1[3]);
            strcat(str, temp_str);
        } else
            for (j = 1; j <= 3; j++)
                strcat(str, bstr);

        ip2 = (parallel) ? phos[2 + p_offset][i + 1] : phos[2 + p_offset][i];  /* strand II */
        if (ip2) {
            ddxyz(mst_org + ioffset3, xyz[ip2], temp);
            for (j = 1; j <= 3; j++)
                P_mst2[j] = dot(temp, &mst_orien[ioffset9 + (j - 1) * 3]);
            if (!parallel) {  /* anti-parallel */
                P_mst2[2] = -P_mst2[2];  /* reverse y */
                P_mst2[3] = -P_mst2[3];  /* reverse z */
            }
            sprintf(temp_str, fmt, P_mst2[1], P_mst2[2], P_mst2[3]);
            strcat(str, temp_str);
        } else
            for (j = 1; j <= 3; j++)
                strcat(str, bstr);

        if (ip1 && ip2)
            avexyz(P_mst1, P_mst2, aveP[i]);

        fprintf(fp, "%s\n", str);
    }
}

/* Following Stephen Harvey:
 *          1  |  Zp - ZpA     chi - chiA |
 *   ABI = --- | ---------- + ----------- |
 *          2  | ZpB - ZpA    chiB - chiA |
 *   where: ZpA = 2.2, ZpB = -0.4
 *          chiA = -157 (203); chiB = -108 (252)
 *   ref: Table 1 of the A-DNA motif paper, JMB2000
*/
static double get_ABI(long idx, double Zp, double **chi_angle)
{
    double ZpA = 2.2, ZpB = -0.4, chiA = 203, chiB = 252;
    double x11, x12, x21, x22, xave, ABI = EMPTY_NUMBER;
    double ZpAB, chiAB, tZp, tchi;

    x11 = chi_angle[1][idx];
    x12 = chi_angle[1][idx + 1];
    x21 = chi_angle[2][idx];
    x22 = chi_angle[2][idx + 1];

    if (x11 > EMPTY_CRITERION && x12 > EMPTY_CRITERION && x21 > EMPTY_CRITERION &&
        x22 > EMPTY_CRITERION) {
        x11 = get_chi360(x11);
        x12 = get_chi360(x12);
        x21 = get_chi360(x21);
        x22 = get_chi360(x22);
        if (in_trans(x11) && in_trans(x12) && in_trans(x21) && in_trans(x22)) {
            xave = (x11 + x12 + x21 + x22) / 4.0;
            ZpAB = ZpB - ZpA;
            chiAB = chiB - chiA;
            tZp = (Zp - ZpA) / ZpAB;
            tchi = (xave - chiA) / chiAB;
            ABI = 0.5 * (tZp + tchi);
        }
    }

    return ABI;
}

/* Calculate and print xyz coordinates of P atoms w.r.t. middle frame
   and middle helical frame. A dinucleotide step is classified as A-,
   B- or TA-like */
void print_PP(long parallel, double **twist_rise, long num_bp, char **bp_seq, long **phos,
              double *mst_orien, double *mst_org, double *mst_orienH, double *mst_orgH,
              double **xyz, long *WC_info, long *bphlx, long abi, long **chi, FILE * fp)
{
    char *bstr = "    --- ", *fmt = "%8.2f%8.2f%8.2f";
    char str[BUF512], **step_info;
    double **aveH, **aveS, **chi_angle = NULL, *ABIval = NULL;
    long i, j, ip, nbpm1, p_offset = 0, ds = 2;
    long *strABT, *idx;

    FILE *fchek;

    nbpm1 = num_bp - 1;

    aveS = dmatrix(1, nbpm1, 1, 3);
    aveH = dmatrix(1, nbpm1, 1, 3);

    fchek = open_file(AUX_FILE, "a");

    /* added August 3, 2004: for O1P & O2P as well */
    p_offset += 2;
    sprintf(str, "xyz coordinates of O1P atoms w.r.t. the middle frame of each dimer\n\n");
    strcat(str, "    step       xI      yI      zI     xII     yII     zII\n");
    print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orien, mst_org, xyz, fchek, str, aveS, p_offset);

    sprintf(str, "xyz coordinates of O1P atoms w.r.t. the middle helix frame"
            " of each dimer\n\n");
    strcat(str, "    step      xIH     yIH     zIH     xIIH    yIIH    zIIH\n");
    print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orienH, mst_orgH, xyz, fchek,
               str, aveH, p_offset);

    p_offset += 2;
    sprintf(str, "xyz coordinates of O2P atoms w.r.t. the middle frame of each dimer\n\n");
    strcat(str, "    step       xI      yI      zI     xII     yII     zII\n");
    print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orien, mst_org, xyz, fchek, str, aveS, p_offset);

    sprintf(str, "xyz coordinates of O2P atoms w.r.t. the middle helix frame"
            " of each dimer\n\n");
    strcat(str, "    step      xIH     yIH     zIH     xIIH    yIIH    zIIH\n");
    print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orienH, mst_orgH, xyz, fchek,
               str, aveH, p_offset);
    /* -------------------------------------------------------------------- */

    p_offset = 0;  /* column index offset */
    sprintf(str, "xyz coordinates of P atoms w.r.t. the middle frame of each dimer\n\n");
    strcat(str, "    step       xI      yI      zI     xII     yII     zII\n");
    print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orien, mst_org, xyz, fchek, str, aveS, p_offset);

    sprintf(str, "xyz coordinates of P atoms w.r.t. the middle helix frame" " of each dimer\n\n");
    strcat(str, "    step      xIH     yIH     zIH     xIIH    yIIH    zIIH\n");
    print_xyzP(parallel, nbpm1, bp_seq, phos, mst_orienH, mst_orgH, xyz, fchek,
               str, aveH, p_offset);

    close_file(fchek);

    if (abi) {
        chi_angle = dmatrix(1, ds, 1, num_bp);
        init_dmatrix(chi_angle, 1, ds, 1, num_bp, EMPTY_NUMBER);
        get_chi_torsions(ds, num_bp, chi, xyz, chi_angle);
        ABIval = dvector(1, num_bp);
        init_dvector(ABIval, 1, num_bp, EMPTY_NUMBER);
    }

    /* classification of each dinucleotide step as A-, B-, TA- or others */
    step_info = cmatrix(1, nbpm1, 0, BUF512);
    strABT = lvector(1, nbpm1);
    idx = lvector(1, num_bp);
    print_sep(fp, '*', 76);
    fprintf(fp, "Classification of each dinucleotide step in"
            " a right-handed nucleic acid\n"
            "structure: A-like; B-like; TA-like; intermediate of A and B, or other cases.\n\n");
    if (abi)
        fprintf(fp, "For definition of the A-B index (ABI), see Waters et al. (2016).\n"
                "``Transitions of Double-Stranded DNA Between the A- and B-Forms.''\n"
                "J. Phys. Chem. B, 120(33), pp8449–8456.\n\n");

    fprintf(fp, "    step       Xp      Yp      Zp     XpH     YpH     ZpH    Form");
    if (abi)
        fprintf(fp, "   ABI");
    fprintf(fp, "\n");

    for (i = 1; i <= nbpm1; i++) {
        sprintf(step_info[i], "%4ld %c%c/%c%c", i, bp_seq[1][i],
                bp_seq[1][i + 1], bp_seq[2][i + 1], bp_seq[2][i]);
        ip = phos[1][i + 1] && ((parallel) ? phos[2][i + 1] : phos[2][i]);
        if (ip && !bphlx[i]) {
            sprintf(str, fmt, aveS[i][1], aveS[i][2], aveS[i][3]);
            strcat(step_info[i], str);
            sprintf(str, fmt, aveH[i][1], aveH[i][2], aveH[i][3]);
            strcat(step_info[i], str);
            if (WC_info[i] && WC_info[i + 1] &&  /* WC geometry */
                dval_in_range(twist_rise[i][1], 10.0, 60.0) &&  /* right-handed */
                dval_in_range(twist_rise[i][2], 2.5, 5.5) &&  /* Rise in range */
                dval_in_range(aveS[i][1], -5.0, -0.5) &&  /* Xp */
                dval_in_range(aveS[i][2], 7.5, 10.0) &&  /* Yp */
                dval_in_range(aveS[i][3], -2.0, 3.5) &&  /* Zp */
                dval_in_range(aveH[i][1], -11.5, 2.5) &&  /* XpH */
                dval_in_range(aveH[i][2], 1.5, 10.0) &&  /* YpH */
                dval_in_range(aveH[i][3], -3.0, 9.0)) {  /* ZpH */
                if (aveS[i][3] >= 1.5)  /* A-form */
                    strABT[i] = 1;
                else if (aveH[i][3] >= 4.0)  /* TA-form */
                    strABT[i] = 3;
                else if (aveS[i][3] <= 0.5 && aveH[i][1] < 0.5)  /* B-form */
                    strABT[i] = 2;  /* aveS[i][3] < 0.5 for C-DNA #47 */
                if (abi)
                    ABIval[i] = get_ABI(i, aveS[i][3], chi_angle);
            }
        } else
            for (j = 1; j <= 7; j++)
                strcat(step_info[i], bstr);
    }
    idx[1] = BUF512;
    idx[num_bp] = BUF512;
    for (i = 2; i <= nbpm1; i++)
        idx[i] = strABT[i] - strABT[i - 1];
    for (i = 1; i <= nbpm1; i++) {
        if (strABT[i] && (num_bp == 2 || !idx[i] || !idx[i + 1])) {
            if (strABT[i] == 1)
                fprintf(fp, "%s     A", step_info[i]);
            else if (strABT[i] == 2)
                fprintf(fp, "%s     B", step_info[i]);
            else
                fprintf(fp, "%s  *TA*", step_info[i]);
        } else
            fprintf(fp, "%s      ", step_info[i]);
        if (abi && ABIval[i] > EMPTY_CRITERION)
            fprintf(fp, " %6.2f", ABIval[i]);
        fprintf(fp, "\n");
    }

    free_dmatrix(aveS, 1, nbpm1, 1, 3);
    free_dmatrix(aveH, 1, nbpm1, 1, 3);
    free_cmatrix(step_info, 1, nbpm1, 0, BUF512);
    free_lvector(strABT, 1, nbpm1);
    free_lvector(idx, 1, num_bp);

    if (abi) {
        free_dmatrix(chi_angle, 1, ds, 1, num_bp);
        free_dvector(ABIval, 1, num_bp);
    }

    /* write P xyz coordinates to "auxiliary.par" */
    print_axyz(num_bp, bp_seq, phos, "P", xyz);
}

/* Classify a structure as right-handed or left handed forms */
void str_classify(double twist_p, double twist_n, long str_type, long parallel,
                  long num_bp, FILE * fp)
{
    long mhlx = 0;

    print_sep(fp, '*', 76);
    fprintf(fp, "Structure classification: \n\n");

    /* added on 2014-dec-02 */
    if (!twist_p && !twist_n && !parallel && str_type == 1) {
        fprintf(fp, "This is a left-handed Z-form structure (as in Z-RNA, 1t4x)\n");
        return;
    }

    if (num_bp == 1) {
        fprintf(fp, "This structure contains only one base-pair\n");
        return;
    }
    if (parallel) {
        fprintf(fp, "This is a parallel duplex structure\n");
        return;
    }
    if (!twist_p && !twist_n)
        return;

    if (str_type >= 10) {
        mhlx = 1;
        fprintf(fp, "This structure contains more than one helical regions\n");
    }

    if (twist_p && twist_n) {
        if (str_type == -2)
            fprintf(fp, "This structure seems to contain a B-Z junction\n");
        if (str_type == -12)
            fprintf(fp, "This structure contains a B-Z junction...\n\n");
        if (str_type < 0) {
            fprintf(fp,
                    "Note that 3DNA determines this junction based *purely* on the geometry of\n"
                    "your input structure. The junction is between a right-handed fragment with\n"
                    "a left-handed one, not necessarily just B-DNA and Z-DNA, even though most\n"
                    "likely that would be the case.\n");
            return;
        }
    }

    str_type %= 10;
    if (str_type == 2) {
        fprintf(fp, "This nucleic acid structure is *unusual*\n");
        return;
    }

    if (!twist_p && twist_n) {
        if (str_type == 1)
            fprintf(fp, "This is a left-handed Z-form structure\n");
        else if (!mhlx)
            fprintf(fp, "This is a left-handed W-form structure\n");
    }

    if (twist_p && !twist_n) {
        fprintf(fp, "This is a right-handed ");
        if (str_type == 1)
            fprintf(fp, "unknown R-form structure\n");
        else
            fprintf(fp, "nucleic acid structure\n");
    }
}

/* Calculate helix radius */
double a_hlxdist(long idx, double **xyz, double *hlx_axis, double *hlx_pos)
{
    double temp, d[4];
    long i;

    if (idx) {
        ddxyz(hlx_pos, xyz[idx], d);
        temp = dot(d, hlx_axis);
        for (i = 1; i <= 3; i++)
            d[i] -= temp * hlx_axis[i];
        return veclen(d);
    } else
        return EMPTY_NUMBER;
}

/* Print the radius from P, O4' & C1' atoms to the local helical axis */
void print_radius(char **bp_seq, long nbpm1, long ich, double **p_radius,
                  double **o4_radius, double **c1_radius, long *bphlx, FILE * fp)
{
    char *bstr = "      ----", *fmt = "%10.2f";
    char str[BUF512];
    long i, j, ik;

    if (!lval_in_range(ich, 1, 3))
        fatal("wrong option [%ld] for printing helix radius\n", ich);

    if (ich == 1) {  /* duplex */
        fprintf(fp, "                        Strand I"
                "                      Strand II\n"
                "     step         P        O4'       C1'" "        P        O4'        C1'\n");
        for (i = 1; i <= nbpm1; i++) {
            sprintf(str, "%4ld %c%c/%c%c", i, bp_seq[1][i],
                    bp_seq[1][i + 1], bp_seq[2][i + 1], bp_seq[2][i]);
            if (bphlx[i])
                for (j = 1; j <= 6; j++)
                    strcat(str, bstr);
            else {
                for (j = 1; j <= 2; j++) {
                    parcat(str, p_radius[j][i], fmt, bstr);
                    parcat(str, o4_radius[j][i], fmt, bstr);
                    parcat(str, c1_radius[j][i], fmt, bstr);
                }
            }
            fprintf(fp, "%s\n", str);
        }

    } else {  /* single strand */
        fprintf(fp, "    step          P        O4'       C1'\n");
        ik = ich - 1;
        for (i = 1; i <= nbpm1; i++) {
            sprintf(str, "%4ld  %c/%c ", i, bp_seq[ik][i], bp_seq[ik][i + 1]);
            parcat(str, p_radius[ik][i], fmt, bstr);
            parcat(str, o4_radius[ik][i], fmt, bstr);
            parcat(str, c1_radius[ik][i], fmt, bstr);
            fprintf(fp, "%s\n", str);
        }
    }
}

/* Get radius from P, O4' and C1' to the helical axis */
void helix_radius(long ds, long num_bp, char **bp_seq, double **orien, double **org,
                  long **phos, long **chi, double **xyz, long *bphlx, FILE * fp)
{
    double temp1, temp2;
    double hx[4], morg[4], o1[44], o2[4], pars[7];
    double **mst, **r1, **r2, **c1_radius, **o4_radius, **p_radius;
    double *bp_orien, *bp_org, **c1BP_radius, **o4BP_radius, **pBP_radius;
    long i, ioffset, j, k, nbpm1, pn;
    FILE *fchek;

    nbpm1 = num_bp - 1;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);

    p_radius = dmatrix(1, ds, 1, nbpm1);
    o4_radius = dmatrix(1, ds, 1, nbpm1);  /* two O4' atom per step */
    c1_radius = dmatrix(1, ds, 1, nbpm1);  /* two C1' atom per step */

    for (i = 1; i <= ds; i++)
        for (j = 1; j <= nbpm1; j++) {
            refs_i_j(j, j + 1, orien[i], org[i], r1, o1, r2, o2);
            helical_par(r1, o1, r2, o2, pars, mst, morg);
            for (k = 1; k <= 3; k++)
                hx[k] = mst[k][3];

            pn = (i == 1) ? j + 1 : j;  /* P index: +1 for I */
            p_radius[i][j] = a_hlxdist(phos[i][pn], xyz, hx, morg);

            ioffset = (j - 1) * 4;

            temp1 = a_hlxdist(chi[i][ioffset + 1], xyz, hx, morg);
            temp2 = a_hlxdist(chi[i][ioffset + 5], xyz, hx, morg);
            o4_radius[i][j] = 0.5 * (temp1 + temp2);

            temp1 = a_hlxdist(chi[i][ioffset + 2], xyz, hx, morg);
            temp2 = a_hlxdist(chi[i][ioffset + 6], xyz, hx, morg);
            c1_radius[i][j] = 0.5 * (temp1 + temp2);
        }

    print_sep(fp, '*', 76);
    fprintf(fp, "Helix radius (radial displacement of P, O4', and C1'"
            " atoms in local helix\n   frame of each dimer)\n\n");

    if (ds == 1)
        print_radius(bp_seq, nbpm1, 2, p_radius, o4_radius, c1_radius, bphlx, fp);
    else {
        bp_org = dvector(1, num_bp * 3);
        bp_orien = dvector(1, num_bp * 9);

        pBP_radius = dmatrix(1, ds, 1, nbpm1);
        o4BP_radius = dmatrix(1, ds, 1, nbpm1);
        c1BP_radius = dmatrix(1, ds, 1, nbpm1);

        for (i = 1; i <= num_bp; i++) {
            refs_right_left(i, orien, org, r1, o1, r2, o2);
            bpstep_par(r1, o1, r2, o2, pars, mst, morg);

            cpxyz(morg, bp_org + (i - 1) * 3);
            mst2orien(bp_orien, (i - 1) * 9, mst);
        }

        for (i = 1; i <= nbpm1; i++) {
            refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
            helical_par(r1, o1, r2, o2, pars, mst, morg);
            for (k = 1; k <= 3; k++)
                hx[k] = mst[k][3];

            ioffset = (i - 1) * 4;

            for (j = 1; j <= ds; j++) {
                pn = (j == 1) ? i + 1 : i;  /* P index: +1 for I */

                pBP_radius[j][i] = a_hlxdist(phos[j][pn], xyz, hx, morg);

                temp1 = a_hlxdist(chi[j][ioffset + 1], xyz, hx, morg);
                temp2 = a_hlxdist(chi[j][ioffset + 5], xyz, hx, morg);
                o4BP_radius[j][i] = 0.5 * (temp1 + temp2);

                temp1 = a_hlxdist(chi[j][ioffset + 2], xyz, hx, morg);
                temp2 = a_hlxdist(chi[j][ioffset + 6], xyz, hx, morg);
                c1BP_radius[j][i] = 0.5 * (temp1 + temp2);
            }
        }
        print_radius(bp_seq, nbpm1, 1, pBP_radius, o4BP_radius, c1BP_radius, bphlx, fp);

        fchek = open_file(AUX_FILE, "a");

        print_sep(fchek, '*', 76);
        fprintf(fchek, "Strand I helix radius\n\n");
        print_radius(bp_seq, nbpm1, 2, p_radius, o4_radius, c1_radius, bphlx, fchek);

        print_sep(fchek, '*', 76);
        fprintf(fchek, "Strand II helix radius\n\n");
        print_radius(bp_seq, nbpm1, 3, p_radius, o4_radius, c1_radius, bphlx, fchek);

        close_file(fchek);

        free_dvector(bp_org, 1, num_bp * 3);
        free_dvector(bp_orien, 1, num_bp * 9);
        free_dmatrix(pBP_radius, 1, ds, 1, nbpm1);
        free_dmatrix(o4BP_radius, 1, ds, 1, nbpm1);
        free_dmatrix(c1BP_radius, 1, ds, 1, nbpm1);
    }

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(p_radius, 1, ds, 1, nbpm1);
    free_dmatrix(o4_radius, 1, ds, 1, nbpm1);
    free_dmatrix(c1_radius, 1, ds, 1, nbpm1);
}

/* Print local base/base-pair helical reference frames */
void print_shlx(char **bp_seq, long nbpm1, long ich, double *shlx_orien,
                double *shlx_org, FILE * fp)
{
    long i, ioffset6, ioffset18, j, k;

    if (!lval_in_range(ich, 1, 3))
        fatal("wrong option [%ld] for printing helix axis\n", ich);

    fprintf(fp, "               Ox      Oy      Oz     Xx    Xy    Xz"
            "   Yx    Yy    Yz    Zx    Zy    Zz\n");

    k = ich - 1;
    for (i = 1; i <= nbpm1; i++) {
        if (ich == 1)
            fprintf(fp, "%4ld %c%c/%c%c", i, bp_seq[1][i],
                    bp_seq[1][i + 1], bp_seq[2][i + 1], bp_seq[2][i]);
        else
            fprintf(fp, "%4ld  %c/%c ", i, bp_seq[k][i], bp_seq[k][i + 1]);

        ioffset6 = (i - 1) * 6;
        ioffset18 = (i - 1) * 18;
        for (j = 1; j <= 3; j++)
            fprintf(fp, "%8.2f", shlx_org[ioffset6 + j]);
        for (j = 1; j <= 9; j++)
            fprintf(fp, "%6.2f", shlx_orien[ioffset18 + j]);
        fprintf(fp, "\n          ");

        ioffset6 += 3;
        ioffset18 += 9;
        for (j = 1; j <= 3; j++)
            fprintf(fp, "%8.2f", shlx_org[ioffset6 + j]);
        for (j = 1; j <= 9; j++)
            fprintf(fp, "%6.2f", shlx_orien[ioffset18 + j]);
        fprintf(fp, "\n");
    }
}

/* Get the local helical axis for bending analysis */
void get_helix_axis(long ds, long num_bp, char **bp_seq, double **orien,
                    double **org, long *bphlx, FILE * fp)
{
    char *bstr = "      ----", *fmt = "%10.2f";
    double hf_twist, hf_rise;
    double hx[4], morg[4], o1[4], o2[4], pars[7];
    double **hlx_org, **hlx_orien, **mst, **r1, **r2, **temp1, **temp2;
    double *bp_hlx_org, *bp_hlx_orien, *bp_org, *bp_orien;
    long i, ioffset, j, k, nbpm1;
    FILE *fchek;

    nbpm1 = num_bp - 1;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);
    temp1 = dmatrix(1, 3, 1, 3);
    temp2 = dmatrix(1, 3, 1, 3);

    hlx_orien = dmatrix(1, ds, 1, nbpm1 * 18);
    hlx_org = dmatrix(1, ds, 1, nbpm1 * 6);

    for (i = 1; i <= ds; i++)
        for (j = 1; j <= nbpm1; j++) {
            refs_i_j(j, j + 1, orien[i], org[i], r1, o1, r2, o2);
            helical_par(r1, o1, r2, o2, pars, mst, morg);
            for (k = 1; k <= 3; k++)
                hx[k] = mst[k][3];

            hf_twist = 0.5 * pars[6];
            hf_rise = 0.5 * pars[3];

            ioffset = (j - 1) * 18;
            arb_rotation(hx, -hf_twist, temp1);
            multi_matrix(temp1, 3, 3, mst, 3, 3, temp2);
            mst2orien(hlx_orien[i], ioffset, temp2);

            ioffset += 9;
            arb_rotation(hx, hf_twist, temp1);
            multi_matrix(temp1, 3, 3, mst, 3, 3, temp2);
            mst2orien(hlx_orien[i], ioffset, temp2);

            ioffset = (j - 1) * 6;
            for (k = 1; k <= 3; k++) {
                hlx_org[i][ioffset + k] = morg[k] - hf_rise * hx[k];
                hlx_org[i][ioffset + k + 3] = morg[k] + hf_rise * hx[k];
            }
        }

    fchek = open_file(AUX_FILE, "a");

    if (ds == 1) {
        print_sep(fp, '*', 76);
        fprintf(fp, "Position (Px, Py, Pz) and local helical axis vector" " (Hx, Hy, Hz)\n\n");
        fprintf(fp, "     step       Px        Py        Pz" "        Hx        Hy        Hz\n");

        for (i = 1; i <= nbpm1; i++) {
            ioffset = (i - 1) * 6;
            avexyz(hlx_org[ds] + ioffset, hlx_org[ds] + ioffset + 3, morg);
            cpxyz(hlx_orien[ds] + (i - 1) * 18 + 6, hx);

            fprintf(fp, "%4ld  %c/%c ", i, bp_seq[ds][i], bp_seq[ds][i + 1]);
            for (j = 1; j <= 3; j++)
                fprintf(fp, fmt, morg[j]);
            for (j = 1; j <= 3; j++)
                fprintf(fp, fmt, hx[j]);
            fprintf(fp, "\n");
        }

        print_sep(fchek, '*', 76);
        fprintf(fchek, "Helix axis\n\n");
        print_shlx(bp_seq, nbpm1, 2, hlx_orien[ds], hlx_org[ds], fchek);

    } else {  /* duplex */
        bp_orien = dvector(1, num_bp * 9);
        bp_org = dvector(1, num_bp * 3);
        bp_hlx_orien = dvector(1, nbpm1 * 18);
        bp_hlx_org = dvector(1, nbpm1 * 6);

        for (i = 1; i <= num_bp; i++) {
            refs_right_left(i, orien, org, r1, o1, r2, o2);
            bpstep_par(r1, o1, r2, o2, pars, mst, morg);

            cpxyz(morg, bp_org + (i - 1) * 3);
            mst2orien(bp_orien, (i - 1) * 9, mst);
        }

        print_sep(fp, '*', 76);
        fprintf(fp, "Position (Px, Py, Pz) and local helical axis vector"
                " (Hx, Hy, Hz)\n         for each dinucleotide step\n\n");
        fprintf(fp, "     step       Px        Py        Pz" "        Hx        Hy        Hz\n");

        for (i = 1; i <= nbpm1; i++) {
            refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
            helical_par(r1, o1, r2, o2, pars, mst, morg);
            for (j = 1; j <= 3; j++)
                hx[j] = mst[j][3];

            hf_twist = 0.5 * pars[6];
            hf_rise = 0.5 * pars[3];

            ioffset = (i - 1) * 18;
            arb_rotation(hx, -hf_twist, temp1);
            multi_matrix(temp1, 3, 3, mst, 3, 3, temp2);
            mst2orien(bp_hlx_orien, ioffset, temp2);

            ioffset += 9;
            arb_rotation(hx, hf_twist, temp1);
            multi_matrix(temp1, 3, 3, mst, 3, 3, temp2);
            mst2orien(bp_hlx_orien, ioffset, temp2);

            ioffset = (i - 1) * 6;
            for (j = 1; j <= 3; j++) {
                bp_hlx_org[ioffset + j] = morg[j] - hf_rise * hx[j];
                bp_hlx_org[ioffset + j + 3] = morg[j] + hf_rise * hx[j];
            }

            fprintf(fp, "%4ld %c%c/%c%c", i, bp_seq[1][i],
                    bp_seq[1][i + 1], bp_seq[2][i + 1], bp_seq[2][i]);
            if (bphlx[i])
                for (j = 1; j <= 6; j++)
                    fprintf(fp, "%s", bstr);
            else {
                for (j = 1; j <= 3; j++)
                    fprintf(fp, fmt, morg[j]);
                for (j = 1; j <= 3; j++)
                    fprintf(fp, fmt, hx[j]);
            }
            fprintf(fp, "\n");
        }

        print_sep(fchek, '*', 88);
        fprintf(fchek, "Base-pair helix axis\n\n");
        print_shlx(bp_seq, nbpm1, 1, bp_hlx_orien, bp_hlx_org, fchek);

        print_sep(fchek, '*', 88);
        fprintf(fchek, "Strand I helix axis\n\n");
        print_shlx(bp_seq, nbpm1, 2, hlx_orien[1], hlx_org[1], fchek);

        print_sep(fchek, '*', 88);
        fprintf(fchek, "Strand II helix axis\n\n");
        print_shlx(bp_seq, nbpm1, 3, hlx_orien[2], hlx_org[2], fchek);

        free_dvector(bp_orien, 1, num_bp * 9);
        free_dvector(bp_org, 1, num_bp * 3);
        free_dvector(bp_hlx_orien, 1, nbpm1 * 18);
        free_dvector(bp_hlx_org, 1, nbpm1 * 6);
    }

    close_file(fchek);

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(temp1, 1, 3, 1, 3);
    free_dmatrix(temp2, 1, 3, 1, 3);
    free_dmatrix(hlx_orien, 1, ds, 1, nbpm1 * 18);
    free_dmatrix(hlx_org, 1, ds, 1, nbpm1 * 6);
}

/* Get the helical axis based on C1* and RN9/YN1 atoms */
void get_axis(long nvec, long **idx, long num, double **xyz, long nb, long *C1b,
              long *C1e, double *std_rise, double *hrise, double *haxis,
              double *hstart, double *hend)
{
    double tb, te, hinge[4], org_xyz[4], t2[3];
    double *g, *drise;
    double **dxy, **dxyT, **rotmat, **vxyz, **xyzH;
    double **dd, **inv_dd;
    long i, j;

    *hrise = EMPTY_NUMBER;
    if (nvec < 3)  /* too few vectors to define a helix */
        return;

    /* find helical axis and rise */
    vxyz = dmatrix(1, nvec, 1, 3);
    drise = dvector(1, nvec);
    for (i = 1; i <= nvec; i++)
        ddxyz(xyz[idx[i][1]], xyz[idx[i][2]], vxyz[i]);
    ls_plane(vxyz, nvec, haxis, hinge, hrise, drise);
    if (*hrise < 0.0) {
        *hrise = -*hrise;
        negate_xyz(haxis);
    }
    *std_rise = std_dvector(drise, nvec);

    /* align haxis to global z-axis */
    rotmat = dmatrix(1, 3, 1, 3);
    xyzH = dmatrix(1, num, 1, 3);
    align2zaxis(num, haxis, rotmat, xyz, xyzH);

    /* locate xy-coordinate the helix passes through */
    dxy = dmatrix(1, nvec, 1, 2);
    dxyT = dmatrix(1, 2, 1, nvec);
    g = dvector(1, nvec);
    for (i = 1; i <= nvec; i++) {
        g[i] = 0.0;
        for (j = 1; j <= 2; j++) {
            tb = xyzH[idx[i][1]][j];
            te = xyzH[idx[i][2]][j];
            dxy[i][j] = 2.0 * (te - tb);
            g[i] += te * te - tb * tb;
        }
    }

    dd = dmatrix(1, 2, 1, 2);
    inv_dd = dmatrix(1, 2, 1, 2);
    multi_vec_matrix(g, nvec, dxy, nvec, 2, t2);
    transpose_matrix(dxy, nvec, 2, dxyT);
    multi_matrix(dxyT, 2, nvec, dxy, nvec, 2, dd);
    dinverse(dd, 2, inv_dd);
    multi_vec_matrix(t2, 2, inv_dd, 2, 2, org_xyz);

    tb = 0.0;  /* ave. z-coordinate of 1st base set C1* atoms */
    j = 0;
    for (i = 1; i <= nb; i++) {
        if (C1b[i]) {
            j++;
            tb += xyzH[C1b[i]][3];
        }
    }
    org_xyz[3] = tb / j;
    multi_vec_matrix(org_xyz, 3, rotmat, 3, 3, hstart);

    te = 0.0;  /* ave. z-coordinate of last base set C1* atoms */
    j = 0;
    for (i = 1; i <= nb; i++) {
        if (C1e[i]) {
            j++;
            te += xyzH[C1e[i]][3];
        }
    }
    org_xyz[3] = te / j;
    multi_vec_matrix(org_xyz, 3, rotmat, 3, 3, hend);

    free_dmatrix(vxyz, 1, nvec, 1, 3);
    free_dvector(drise, 1, nvec);
    free_dmatrix(rotmat, 1, 3, 1, 3);
    free_dmatrix(xyzH, 1, num, 1, 3);
    free_dmatrix(dxy, 1, nvec, 1, 2);
    free_dmatrix(dxyT, 1, 2, 1, nvec);
    free_dvector(g, 1, nvec);
    free_dmatrix(dd, 1, 2, 1, 2);
    free_dmatrix(inv_dd, 1, 2, 1, 2);
}

/* Print out P, O4' and C1' atom radius defined cylinder in Raster3D format */
void print_poc_r3d(double *rave, double *hstart, double *hend)
{
    char label_style[BUF512];
    double width3[4], hb_col[5], **atom_col, **base_col;
    long i;
    FILE *fpr3d;

    atom_col = dmatrix(0, NATOMCOL, 1, 3);
    base_col = dmatrix(0, NBASECOL, 1, 3);
    get_r3dpars(base_col, hb_col, width3, atom_col, label_style);

    fpr3d = open_file(POC_FILE, "w");
    fprintf(fpr3d, "###\n### Linear global helical axis if not strongly curved\n");

    rave[3] = width3[1];
    for (i = 0; i < 4; i++)  /* 0 for P; 1 for O, 2 for C, and 3 for thinner line */
        r3d_rod((i != 3) ? 99L : 5L, hstart, hend, rave[i], hb_col, fpr3d);
    close_file(fpr3d);

    free_dmatrix(atom_col, 0, NATOMCOL, 1, 3);
    free_dmatrix(base_col, 0, NBASECOL, 1, 3);
}

/* Get the radius of P, O4' and C1' atoms from the global helical axis */
static void get_poc_radius(long ds, long num_bp, long **phos, long **chi, double **xyz,
                           double *haxis, double *hstart, double *hend, double *rave, double *rstd)
{
    long i, ik, ioffset, j, k, num_ple, poc_num[3];
    double **poc_radius;

    for (i = 0; i <= 2; i++) {
        poc_num[i] = 0;
        rave[i] = rstd[i] = 0.0;
    }

    num_ple = 2 * num_bp;  /* maximum possible number */
    poc_radius = dmatrix(0, 2, 1, num_ple);
    for (i = 1; i <= ds; i++)
        for (j = 1; j <= num_bp; j++) {
            k = 0;  /* P */
            ik = phos[i][j];
            if (ik)
                poc_radius[k][++poc_num[k]] = a_hlxdist(ik, xyz, haxis, hstart);
            ioffset = (j - 1) * 4;
            k = 1;  /* O4 */
            ik = chi[i][ioffset + 1];
            if (ik)
                poc_radius[k][++poc_num[k]] = a_hlxdist(ik, xyz, haxis, hstart);
            k = 2;  /* C1 */
            ik = chi[i][ioffset + 2];
            if (ik)
                poc_radius[k][++poc_num[k]] = a_hlxdist(ik, xyz, haxis, hstart);
        }

    for (i = 0; i <= 2; i++)
        if (poc_num[i]) {
            rave[i] = ave_dvector(poc_radius[i], poc_num[i]);
            rstd[i] = std_dvector(poc_radius[i], poc_num[i]);
        }
    print_poc_r3d(rave, hstart, hend);

    free_dmatrix(poc_radius, 0, 2, 1, num_ple);
}

/* Global helical parameters based on C1'-C1' vector for a duplex */
static void global_c1_c1_par(long num_bp, char **bp_seq, long **chi, double **xyz,
                             double *haxis, double *hstart, FILE * fp)
{
    char *bstr = "      --- ", *fmt = "%10.2f", **c1_hpar;
    double dtmp, dd[4], org_xyz[4], **dc1, **mc1;
    long i, ia, ib, ioffset, j;

    fprintf(fp, "\nGlobal parameters based on C1'-C1' vectors:\n\n");
    fprintf(fp,
            "disp.: displacement of the middle C1'-C1' point from the helix\n"
            "angle: inclination between C1'-C1' vector and helix (subtracted from 90)\n"
            "twist: helical twist angle between consecutive C1'-C1' vectors\n"
            "rise:  helical rise by projection of the vector connecting consecutive\n"
            "       C1'-C1' middle points onto the helical axis\n\n");

    c1_hpar = cmatrix(1, num_bp, 0, 50);
    dc1 = dmatrix(1, num_bp, 1, 3);
    mc1 = dmatrix(1, num_bp, 1, 3);
    for (i = 1; i <= num_bp; i++) {  /* displacement and angle */
        sprintf(c1_hpar[i], "%4ld %c%c%c", i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
        ioffset = (i - 1) * 4;
        ia = chi[1][ioffset + 2];
        ib = chi[2][ioffset + 2];
        if (ia && ib) {
            ddxyz(xyz[ib], xyz[ia], dc1[i]);
            avexyz(xyz[ia], xyz[ib], mc1[i]);
            ddxyz(hstart, mc1[i], dd);
            dtmp = dot(dd, haxis);
            for (j = 1; j <= 3; j++)
                org_xyz[j] = dd[j] - dtmp * haxis[j];
            parcat(c1_hpar[i], veclen(org_xyz), fmt, bstr);
            dtmp = 90.0 - magang(dc1[i], haxis);
            parcat(c1_hpar[i], dtmp, fmt, bstr);
        } else {
            parcat(c1_hpar[i], EMPTY_NUMBER, fmt, bstr);
            parcat(c1_hpar[i], EMPTY_NUMBER, fmt, bstr);
        }
    }

    for (i = 1; i <= num_bp - 1; i++) {  /* twist and rise */
        if (strstr(c1_hpar[i], bstr) == NULL && strstr(c1_hpar[i + 1], bstr) == NULL) {
            dtmp = vec_ang(dc1[i], dc1[i + 1], haxis);
            parcat(c1_hpar[i], dtmp, fmt, bstr);
            ddxyz(mc1[i], mc1[i + 1], dd);
            parcat(c1_hpar[i], dot(dd, haxis), fmt, bstr);
        } else {
            parcat(c1_hpar[i], EMPTY_NUMBER, fmt, bstr);
            parcat(c1_hpar[i], EMPTY_NUMBER, fmt, bstr);
        }
    }
    parcat(c1_hpar[num_bp], EMPTY_NUMBER, fmt, bstr);
    parcat(c1_hpar[num_bp], EMPTY_NUMBER, fmt, bstr);

    fprintf(fp, "     bp       disp.    angle     twist      rise\n");
    for (i = 1; i <= num_bp; i++)
        fprintf(fp, "%s\n", c1_hpar[i]);

    free_cmatrix(c1_hpar, 1, num_bp, 0, 50);
    free_dmatrix(dc1, 1, num_bp, 1, 3);
    free_dmatrix(mc1, 1, num_bp, 1, 3);
}

/* Structural analysis from a global prospective of the backbone:
   measuring curvature, and bending angle */
void global_analysis(long ds, long num_bp, long num, char **bp_seq, long **chi,
                     long **phos, double **xyz, FILE * fp)
{
    double hrise, std_rise, haxis[4], hstart[4], hend[4], rave[4], rstd[3];
    long nvec, C1b[MBASES], C1e[MBASES], **idx;

    idx = lmatrix(1, 4 * num_bp, 1, 2);  /* beginning & end index */
    get_CNidx(ds, num_bp, chi, idx, &nvec, C1b, C1e);
    get_axis(nvec, idx, num, xyz, ds, C1b, C1e, &std_rise, &hrise, haxis, hstart, hend);
    free_lmatrix(idx, 1, 4 * num_bp, 1, 2);

    if (hrise < EMPTY_CRITERION)
        return;  /* no helix defined */

    print_sep(fp, '*', 76);
    fprintf(fp, "Global linear helical axis defined by equivalent C1'"
            " and RN9/YN1 atom pairs\n");
    fprintf(fp, "Deviation from regular linear helix: %.2f(%.2f)\n", hrise, std_rise);

    if (std_rise > Gvars.misc_pars.std_curved)
        return;

    get_poc_radius(ds, num_bp, phos, chi, xyz, haxis, hstart, hend, rave, rstd);

    fprintf(fp, "Helix:  %9.4f %9.4f %9.4f\n", haxis[1], haxis[2], haxis[3]);
    fprintf(fp, "HETATM 9998  XS    X X 999    %8.3f%8.3f%8.3f\n",
            hstart[1], hstart[2], hstart[3]);
    fprintf(fp, "HETATM 9999  XE    X X 999    %8.3f%8.3f%8.3f\n", hend[1], hend[2], hend[3]);
    fprintf(fp, "Average and standard deviation of helix radius:\n");
    fprintf(fp, "      P: %.2f(%.2f), O4': %.2f(%.2f),  C1': %.2f(%.2f)\n",
            rave[0], rstd[0], rave[1], rstd[1], rave[2], rstd[2]);

    if (ds == 1)
        return;

    global_c1_c1_par(num_bp, bp_seq, chi, xyz, haxis, hstart, fp);
}

/* Check if two nucleotides along single strand are reasonably stacked [oct-12-2007] */
static long with_ss_overlap(long istep, long r1, long r2, double **orien,
                            long **ring_atom, double **xyz)
{
    long i, ioffset, j, n1, n2;
    double sdist = XBIG, oave[4] = { EMPTY_NUMBER, 0.0, 0.0, 0.0 };
    double **oxyz1, **oxyz2;

    ioffset = (istep - 1) * 9;

    if (z1_z2_angle_in_0_to_90(&orien[1][ioffset + 6], &orien[1][ioffset + 15]) > 65.0)
        return 0;

    oxyz1 = dmatrix(1, 9, 1, 3);
    oxyz2 = dmatrix(1, 9, 1, 3);

    n1 = ratom_xyz(ring_atom[r1], 0, xyz, oave, oxyz1);  /* including exocyclic atoms */
    n2 = ratom_xyz(ring_atom[r2], 0, xyz, oave, oxyz2);

    for (i = 1; i <= n1; i++)
        for (j = 1; j <= n2; j++)
            sdist = dval_min(sdist, p1p2_dist(oxyz1[i], oxyz2[j]));

    free_dmatrix(oxyz1, 1, 9, 1, 3);
    free_dmatrix(oxyz2, 1, 9, 1, 3);

    if (sdist > 4.5)
        return 0;

    return 1;
}

/* added 2012-12-29 to correct the asymmetry in 3ubt, model A (chains E/F)
 * step #3: GC/GC is intercalated along strand E, not F. Check both chains
 * to ensure consistent results for base-stacking areas. */
static long with_stacking(long ds, long istep, long **pair_num, double **orien,
                          long **ring_atom, double **xyz)
{
    long i1, i2, k1, k2;

    i1 = pair_num[1][istep];
    i2 = pair_num[1][istep + 1];

    if (ds == 1)
        return with_ss_overlap(istep, i1, i2, orien, ring_atom, xyz);

    k1 = pair_num[ds][istep];
    k2 = pair_num[ds][istep + 1];

    return with_ss_overlap(istep, i1, i2, orien, ring_atom, xyz) &&
        with_ss_overlap(istep, k2, k1, orien, ring_atom, xyz);
}

/* Get overlap area between neighbor bases (ring + first order of shell) */
void base_overlap(long ds, long num_bp, long num, long num_residue, long **pair_num,
                  long *bRY, char **bp_seq, long **seidx, char **AtomName, double **xyz,
                  long *idx, double **orien, double **org, FILE * fp)
{
    char *fmt = " %5.2f(%5.2f)";
    long i, j, k, r1, r2;
    long **ring_atom;
    double oave[4], zave[4], *olarea;

    /* 1-9 ring atom index, 10 # of ring atoms, 11-19 first level */
    ring_atom = lmatrix(1, num_residue, 1, 19);
    ring_oidx(num, num_residue, bRY, seidx, AtomName, xyz, idx, ring_atom);

    /* i1-i2(1), i1-j2(2), j1-i2(3), j1-j2(4) ==> sum(5): ring + 1st level */
    olarea = dvector(1, 10);

    print_sep(fp, '*', 76);
    fprintf(fp, "Overlap area in Angstrom^2 between polygons defined by atoms on"
            " successive\nbases. Polygons projected in the mean plane of the designed"
            " base-pair step.\n\n");
    fprintf(fp, "Values in parentheses measure the overlap of base ring atoms only."
            " Those\noutside parentheses include exocyclic atoms on the ring. Intra-"
            " and\ninter-strand overlap is designated according to the following" " diagram:\n\n");
    if (ds == 2) {
        fprintf(fp,
                "                    i2  3'      5' j2\n"
                "                       /|\\      |\n"
                "                        |       |\n"
                "               Strand I |       | II\n"
                "                        |       |\n"
                "                        |       |\n"
                "                        |      \\|/\n"
                "                    i1  5'      3' j1\n\n");
        fprintf(fp, "     step      i1-i2        i1-j2        j1-i2        j1-j2" "        sum\n");
    } else {

        fprintf(fp,
                "                    i2  3'\n"
                "                       /|\\\n"
                "                        |\n"
                "               Strand I |\n"
                "                        |\n"
                "                        |\n"
                "                        |\n" "                    i1  5'\n\n");
        fprintf(fp, "     step      i1-i2\n");
    }

    k = (ds == 1) ? 1 : 4;
    for (i = 1; i < num_bp; i++) {
        init_dvector(olarea, 1, 10, 0.0);

        get_zoave(i, ds, orien, org, oave, zave);

        r1 = pair_num[1][i];  /* i1 */
        r2 = pair_num[1][i + 1];  /* i2 */

        if (with_stacking(ds, i, pair_num, orien, ring_atom, xyz)) {
            olarea[1] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 1);
            olarea[6] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 0);

            if (ds == 2) {
                r2 = pair_num[ds][i + 1];  /* j2 */
                olarea[2] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 1);
                olarea[7] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 0);
                r1 = pair_num[ds][i];  /* j1 */
                olarea[4] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 1);
                olarea[9] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 0);
                r2 = pair_num[1][i + 1];  /* i2 */
                olarea[3] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 1);
                olarea[8] = get_oarea(r1, r2, ring_atom, oave, zave, xyz, 0);
            }

            for (j = 1; j <= k; j++) {
                olarea[5] += olarea[j];
                olarea[10] += olarea[5 + j];
            }
        }

        if (ds == 1)
            fprintf(fp, "%4ld  %c/%c ", i, bp_seq[ds][i], bp_seq[ds][i + 1]);
        else
            fprintf(fp, "%4ld %c%c/%c%c", i, bp_seq[1][i], bp_seq[1][i + 1],
                    bp_seq[ds][i + 1], bp_seq[ds][i]);
        for (j = 1; j <= k; j++)
            fprintf(fp, fmt, olarea[5 + j], olarea[j]);
        if (ds == 2)
            fprintf(fp, fmt, olarea[10], olarea[5]);
        fprintf(fp, "\n");
    }

    free_lmatrix(ring_atom, 1, num_residue, 1, 19);
    free_dvector(olarea, 1, 10);
}

/* Get xyz coordinates for base ring + 1st shell, offset by "oave" */
long ratom_xyz(long *ratom_list, long only_ring, double **xyz, double *oave, double **oxyz)
{
    long i, k, n;

    n = ratom_list[10];
    for (i = 1; i <= n; i++) {
        k = (only_ring) ? ratom_list[i] : ratom_list[10 + i];
        ddxyz(oave, xyz[k], oxyz[i]);
    }

    return n;
}

/* Get average origin and z-axis for program "analyze" */
void get_zoave(long istep, long ds, double **orien, double **org, double *oave, double *zave)
{
    double d;
    long ioffset3, ioffset9, j;

    ioffset3 = (istep - 1) * 3;
    ioffset9 = (istep - 1) * 9;

    d = dot(&orien[1][ioffset9 + 6], &orien[1][ioffset9 + 15]);  /* i1 vs i2 */
    avexyz(org[1] + ioffset3, org[1] + ioffset3 + 3, oave);
    (d > 0.0) ? sumxyz(orien[1] + ioffset9 + 15, orien[1] + ioffset9 + 6, zave) :
        ddxyz(orien[1] + ioffset9 + 15, orien[1] + ioffset9 + 6, zave);

    if (ds == 2) {
        for (j = 1; j <= 3; j++) {
            oave[j] += 0.5 * (org[ds][ioffset3 + j] + org[ds][ioffset3 + 3 + j]);
            oave[j] *= 0.5;
        }
        d = dot(&orien[1][ioffset9 + 6], &orien[ds][ioffset9 + 6]);  /* i1 vs j1 */
        for (j = 1; j <= 3; j++)
            zave[j] += (d > 0.0) ? orien[ds][ioffset9 + 6 + j] : -orien[ds][ioffset9 + 6 + j];
        d = dot(&orien[1][ioffset9 + 6], &orien[ds][ioffset9 + 15]);  /* i1 vs j2 */
        for (j = 1; j <= 3; j++)
            zave[j] += (d > 0.0) ? orien[ds][ioffset9 + 15 + j] : -orien[ds][ioffset9 + 15 + j];
    }

    vec_norm(zave);
}

/* Get average origin and z-axis for a base-pair */
void get_bp_zoave(long ia, long ib, double **orien, double **org, double *oave, double *zave)
{
    double d;

    d = dot(&orien[ia][6], &orien[ib][6]);
    avexyz(org[ia], org[ib], oave);
    (d > 0.0) ? sumxyz(orien[ib] + 6, orien[ia] + 6, zave) :
        ddxyz(orien[ib] + 6, orien[ia] + 6, zave);
    vec_norm(zave);
}

/* Get ring atom index and linkage information in base residues: exclude Hs */
void ring_oidx(long num, long num_residue, long *RY, long **seidx, char **AtomName,
               double **xyz, long *idx, long **ring_atom)
{
    long i, num_ring;
    long **connect;

    connect = lmatrix(1, num, 1, 7);  /* overall connection */
    get_bonds(num, AtomName, xyz, num_residue, RY, seidx, connect);
    all_bring_atoms(num_residue, RY, seidx, AtomName, &num_ring, ring_atom);

    for (i = 1; i <= num_residue; i++)  /* for each residue */
        if (ring_atom[i][10] > 0)
            get_cntatom(ring_atom[i], connect, idx);

    free_lmatrix(connect, 1, num, 1, 7);
}

/* Get ring atom connections */
void get_cntatom(long *ringlist, long **connect, long *idx)
{
    long i, ic = 0, id, ix, j, ra_num;

    ra_num = ringlist[10];  /* # of ring atoms */
    for (i = 1; i <= ra_num; i++) {
        id = ringlist[i];
        ix = 0;
        for (j = 1; j <= connect[id][7]; j++) {
            ic = connect[id][j];
            if (idx[ic] == 3)
                continue;  /* excluding H atom */
            if (!lval_in_set(ic, 1, ra_num, ringlist)) {  /* not a ring atom */
                ix = 1;
                break;
            }
        }
        ringlist[10 + i] = (ix) ? ic : id;
    }
}

/* ============ calculating polygon overlap area ============ */
static void pia_fit(double minx, double miny, const double mid, double sclx, double scly,
                    point * x, long cx, vertex * ix, long fudge)
{
    long c, t;
    point t1, t2;

    c = cx;
    while (c--) {
        t = (x[c].x - minx) * sclx - mid;
        ix[c].ip.x = (t & ~7) | fudge | (c & 1);
        t = (x[c].y - miny) * scly - mid;
        ix[c].ip.y = (t & ~7) | fudge;
    }

    ix[0].ip.y += cx & 1;
    ix[cx] = ix[0];

    c = cx;
    while (c--) {
        t1.x = ix[c].ip.x;
        t1.y = ix[c + 1].ip.x;
        t2.x = ix[c + 1].ip.x;
        t2.y = ix[c].ip.x;
        ix[c].rx = (ix[c].ip.x < ix[c + 1].ip.x) ? t1 : t2;

        t1.x = ix[c].ip.y;
        t1.y = ix[c + 1].ip.y;
        t2.x = ix[c + 1].ip.y;
        t2.y = ix[c].ip.y;
        ix[c].ry = (ix[c].ip.y < ix[c + 1].ip.y) ? t1 : t2;
        ix[c].inside = 0;
    }
}

static double pia_area(point a, point p, point q)
{
    return p.x * q.y - p.y * q.x + a.x * (p.y - q.y) + a.y * (q.x - p.x);
}

static void pia_cntrib(double *s, point f, point t, long w)
{
    *s += w * (t.x - f.x) * (t.y + f.y) * 0.5;
}

static long pia_ovl(point p, point q)
{
    return p.x < q.y && q.x < p.y;
}

static void pia_cross(double *out_s, vertex * a, vertex * b, vertex * c, vertex * d,
                      double a1, double a2, double a3, double a4)
{
    double r1, r2;
    point dp;

    r1 = a1 / (a1 + a2);
    r2 = a3 / (a3 + a4);

    dp.x = a->ip.x + r1 * (b->ip.x - a->ip.x);
    dp.y = a->ip.y + r1 * (b->ip.y - a->ip.y);
    pia_cntrib(out_s, dp, b->ip, 1);

    dp.x = c->ip.x + r2 * (d->ip.x - c->ip.x);
    dp.y = c->ip.y + r2 * (d->ip.y - c->ip.y);
    pia_cntrib(out_s, d->ip, dp, 1);

    ++a->inside;
    --c->inside;
}

static void pia_inness(double *out_s, vertex * P, long cP, vertex * Q, long cQ)
{
    long j, sgn, s = 0, c = cQ;
    point p = P[0].ip;

    while (c--)
        if (Q[c].rx.x < p.x && p.x < Q[c].rx.y) {
            sgn = 0 < pia_area(p, Q[c].ip, Q[c + 1].ip);
            s += sgn != (Q[c].ip.x < Q[c + 1].ip.x) ? 0 : (sgn ? -1 : 1);
        }
    for (j = 0; j < cP; ++j) {
        if (s)
            pia_cntrib(out_s, P[j].ip, P[j + 1].ip, s);
        s += P[j].inside;
    }
}

/* Area of intersection between two polygons */
static double pia_inter(point * a, long na, point * b, long nb)
{
    double minx, miny, maxx, maxy, ascale, sclx, scly;
    double out_s = 0.0, a1, a2, a3, a4;
    const double gamut = 5.0e8, mid = 0.5 * gamut;
    long j, k, o;
    vertex ipa[MNPOLY + 1], ipb[MNPOLY + 1];  /* ISO C89 forbids variable-size array */

    if (na < 3 || nb < 3)
        return 0;

    minx = miny = XBIG;
    maxx = maxy = -XBIG;
    for (j = 0; j < na; j++) {
        if (minx > a[j].x)
            minx = a[j].x;
        if (miny > a[j].y)
            miny = a[j].y;
        if (maxx < a[j].x)
            maxx = a[j].x;
        if (maxy < a[j].y)
            maxy = a[j].y;
    }
    for (j = 0; j < nb; j++) {
        if (minx > b[j].x)
            minx = b[j].x;
        if (miny > b[j].y)
            miny = b[j].y;
        if (maxx < b[j].x)
            maxx = b[j].x;
        if (maxy < b[j].y)
            maxy = b[j].y;
    }

    sclx = gamut / (maxx - minx);
    scly = gamut / (maxy - miny);
    ascale = sclx * scly;

    pia_fit(minx, miny, mid, sclx, scly, a, na, ipa, 0);
    pia_fit(minx, miny, mid, sclx, scly, b, nb, ipb, 2);

    for (j = 0; j < na; ++j)
        for (k = 0; k < nb; ++k)
            if (pia_ovl(ipa[j].rx, ipb[k].rx)
                && pia_ovl(ipa[j].ry, ipb[k].ry)) {
                a1 = -pia_area(ipa[j].ip, ipb[k].ip, ipb[k + 1].ip);
                a2 = pia_area(ipa[j + 1].ip, ipb[k].ip, ipb[k + 1].ip);
                o = a1 < 0;
                if (o == (a2 < 0)) {
                    a3 = pia_area(ipb[k].ip, ipa[j].ip, ipa[j + 1].ip);
                    a4 = -pia_area(ipb[k + 1].ip, ipa[j].ip, ipa[j + 1].ip);
                    if ((a3 < 0) == (a4 < 0)) {
                        if (o)
                            pia_cross(&out_s, &ipa[j], &ipa[j + 1],
                                      &ipb[k], &ipb[k + 1], a1, a2, a3, a4);
                        else
                            pia_cross(&out_s, &ipb[k], &ipb[k + 1],
                                      &ipa[j], &ipa[j + 1], a3, a4, a1, a2);
                    }
                }
            }

    pia_inness(&out_s, ipa, na, ipb, nb);
    pia_inness(&out_s, ipb, nb, ipa, na);

    return fabs(out_s) / ascale;
}

/* ============ calculating polygon overlap area ============ */

/* Get overlap area between base residues r1 and r2 */
double get_oarea(long r1, long r2, long **ring_atom, double *oave, double *zave,
                 double **xyz, long only_ring)
{
    long i, j, n1, n2;
    double **oxyz1, **oxyz2, **oxyz1Z, **oxyz2Z, **rotmat;
    point a[MNPOLY], b[MNPOLY];

    oxyz1 = dmatrix(1, 9, 1, 3);
    oxyz2 = dmatrix(1, 9, 1, 3);
    oxyz1Z = dmatrix(1, 9, 1, 3);
    oxyz2Z = dmatrix(1, 9, 1, 3);
    rotmat = dmatrix(1, 3, 1, 3);  /* no use here */

    n1 = ratom_xyz(ring_atom[r1], only_ring, xyz, oave, oxyz1);
    n2 = ratom_xyz(ring_atom[r2], only_ring, xyz, oave, oxyz2);

    align2zaxis(n1, zave, rotmat, oxyz1, oxyz1Z);
    align2zaxis(n2, zave, rotmat, oxyz2, oxyz2Z);

    /* change xy1 & xy2 to an array of point */
    for (i = 1; i <= n1; i++) {
        j = i - 1;
        a[j].x = oxyz1Z[i][1];
        a[j].y = oxyz1Z[i][2];
    }
    for (i = 1; i <= n2; i++) {
        j = i - 1;
        b[j].x = oxyz2Z[i][1];
        b[j].y = oxyz2Z[i][2];
    }

    free_dmatrix(oxyz1, 1, 9, 1, 3);
    free_dmatrix(oxyz2, 1, 9, 1, 3);
    free_dmatrix(oxyz1Z, 1, 9, 1, 3);
    free_dmatrix(oxyz2Z, 1, 9, 1, 3);
    free_dmatrix(rotmat, 1, 3, 1, 3);

    return pia_inter(a, n1, b, n2);
}

void verify_oarea(void)
{
    long na, nb, nc;
    point a[] = {
        {
         0, 0},
        {
         10, 0},
        {
         10, 10},
        {
         0, 10}
    };
    point b[] = {
        {
         5, 0},
        {
         10, 0},
        {
         10, 5},
        {
         5, 5}
    };
    point c[] = {
        {
         -3, -2},
        {
         -1, 4},
        {
         6, 1},
        {
         3, 10},
        {
         -4, 9},
    };

    na = NELEMS(a);
    nb = NELEMS(b);
    nc = NELEMS(c);

    printf("Area(%ld): %g\n", na, pia_inter(a, na, a, na));
    printf("Area(%ld): %g\n", nb, pia_inter(b, nb, b, nb));
    printf("Area(%ld): %g\n", nc, pia_inter(c, nc, c, nc));

    printf("Overlap area: %g\n", pia_inter(a, na, b, nb));
    printf("Overlap area: %g\n", pia_inter(a, na, c, nc));
    printf("Overlap area: %g\n", pia_inter(b, nb, c, nc));
}

/* Get non-hydrogen base atom index in a residue based on atom name */
void cehs_base_atoms(char **AtomName, long ib, long ie, long *num_batom, long *batom)
{
    long i;

    *num_batom = 0;

    for (i = ib; i <= ie; i++)
        if (is_baseatom(AtomName[i])) {
            if (++*num_batom > NUM_BASE_ATOMS)
                fatal("too many base atoms in a residue\n");
            batom[*num_batom] = i;
        }
}

/* Calculate the six local CEHS base-pair parameters.
   Propeller is applied first followed by buckle-opening */
void cehs_bppar(double **rot1, double *org1, double **rot2, double *org2, double *pars,
                double **mst_orien, double *mst_org)
{
    double buckleopening, phi;
    double hinge[4], t1[4], t2[4], xm[4], ym[4], zm[4];
    double **paraII, **paraI, **temp;
    long i, j;

    for (i = 1; i <= 3; i++) {
        t1[i] = rot1[i][2];  /* y1 */
        t2[i] = rot2[i][2];  /* y2 */
    }

    cross(t1, t2, hinge);
    buckleopening = magang(t1, t2);

    temp = dmatrix(1, 3, 1, 3);
    paraII = dmatrix(1, 3, 1, 3);
    paraI = dmatrix(1, 3, 1, 3);

    arb_rotation(hinge, -0.5 * buckleopening, temp);
    multi_matrix(temp, 3, 3, rot2, 3, 3, paraI);
    arb_rotation(hinge, 0.5 * buckleopening, temp);
    multi_matrix(temp, 3, 3, rot1, 3, 3, paraII);

    for (i = 1; i <= 3; i++) {
        ym[i] = paraI[i][2];  /* also paraII[i][2] */
        t1[i] = paraII[i][1];  /* x1 */
        t2[i] = paraI[i][1];  /* x2 */
    }

    /* twist is the angle between the two y- or x-axes */
    pars[5] = vec_ang(t1, t2, ym);

    sumxyz(t1, t2, xm);
    vec_norm(xm);

    cross(xm, ym, zm);

    avexyz(org1, org2, mst_org);
    ddxyz(org1, org2, t1);
    x_y_z_2_mtx(xm, ym, zm, mst_orien);

    /* get the xyz displacement parameters */
    for (i = 1; i <= 3; i++) {
        pars[i] = 0.0;
        for (j = 1; j <= 3; j++)
            pars[i] += t1[j] * mst_orien[j][i];
    }

    /* phi angle is defined by hinge and xm */
    phi = deg2rad(vec_ang(hinge, xm, ym));

    /* get buckle and opening angles */
    pars[4] = buckleopening * cos(phi);
    pars[6] = buckleopening * sin(phi);

    free_dmatrix(temp, 1, 3, 1, 3);
    free_dmatrix(paraII, 1, 3, 1, 3);
    free_dmatrix(paraI, 1, 3, 1, 3);
}

/* 3DNA uses the CEHS/SCHNAaP scheme but based on a newly recommended reference
   frame which is base-centered as in Curves, RNA and CompDNA. On the other hand,
   CEHS, NUPARM and FreeHelix/NewHelix use RC8 and YC6 atom to define base-pair
   origin and long axis. They represent two classes of reference frames. This
   C6/C8 style parameters can be compared to the 3DNA ones to show how the same
   mathematics could give different results depending on the choice of reference
   frame. */
void cehs_pars(long num_bp, long istart, long istep, long **pair_num, char **bp_seq,
               long **seidx, long **c6_c8, long *RY, char **AtomName, char **ResName,
               char *ChainID, long *ResSeq, char **Miscs, double **xyz, double *bp_orien,
               double *bp_org, long bz, FILE * fp)
{
    char idmsg[BUF512], pn1[5], pn2[5];
    char **step_str;

    double odist;
    double adist[2 * NUM_BASE_ATOMS];
    double dorg[4], o1[4], o2[4], ppos[4], y[4], z[4], mfoi[4];
    double **bxyz, **org, **orien, **mfi, **r1, **r2;
    double **bp_par, **step_par;

    long ds = 2, num_step = 0, bz_junction = 0, z_step = 0;
    long i, ib, ie, ik, ioffset, j, k, rnum, rnum2;
    long ib2, ie2, nbpm1, num_batom, num_batom2, p1, p2, p3;
    long batom[NUM_BASE_ATOMS], batom2[NUM_BASE_ATOMS];
    long *bphlx;  /* always count continously */

    nbpm1 = num_bp - 1;
    bxyz = dmatrix(1, 2 * NUM_BASE_ATOMS, 1, 3);

    /* y-, z- and origin for each base */
    orien = dmatrix(1, ds, 1, num_bp * 9);
    org = dmatrix(1, ds, 1, num_bp * 3);

    for (i = 1; i <= ds; i++)
        for (j = 1; j <= num_bp; j++) {
            rnum = pair_num[i][j];
            ib = seidx[rnum][1];
            ie = seidx[rnum][2];

            get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);

            if (RY[rnum] == 1) {
                strcpy(pn1, " C4 ");
                strcpy(pn2, " N1 ");
            } else if (RY[rnum] == 0) {
                strcpy(pn1, " C6 ");
                strcpy(pn2, " N3 ");
            } else
                fatal("Non-base residue: %s\n", ResName[ib]);

            p1 = find_1st_atom(pn1, AtomName, ib, ie, idmsg);
            p2 = find_1st_atom(pn2, AtomName, ib, ie, idmsg);
            if (RY[rnum] == 0 && (bp_seq[i][j] == 'P' || bp_seq[i][j] == 'p'))
                p3 = find_1st_atom(" C4 ", AtomName, ib, ie, idmsg);
            else
                p3 = find_1st_atom(" C2 ", AtomName, ib, ie, idmsg);
            if (!p1 || !p2 || !p3 || !c6_c8[i][j])
                fatal("missing required atom(s) in %s\n", idmsg);
            avexyz(xyz[p1], xyz[p2], org[i] + (j - 1) * 3);
            ddxyz(xyz[p2], xyz[p1], y);  /* y-axis */
            vec_norm(y);
            if (i == 2)  /* strand II */
                negate_xyz(y);
            ioffset = (j - 1) * 9;
            cpxyz(y, orien[i] + ioffset + 3);

            /* base atoms. use ls-plane normal as z-axis */
            cehs_base_atoms(AtomName, ib, ie, &num_batom, batom);
            for (k = 1; k <= num_batom; k++)
                cpxyz(xyz[batom[k]], bxyz[k]);
            ls_plane(bxyz, num_batom, z, ppos, &odist, adist);

            /* define the direction of z */
            ddxyz(xyz[p2], xyz[p1], o1);
            ddxyz(xyz[p2], xyz[p3], o2);
            cross(o1, o2, ppos);  /* along 5'-->3' of strand I */
            bp_seq[0][j] = (i == 2 && dot(ppos, &orien[1][ioffset + 6]) < 0.0) ? '-' : '+';

            if ((i == 1 && dot(z, ppos) < 0) || (i == 2 && dot(z, ppos) > 0.0))
                negate_xyz(z);

            vec_orth(z, y);  /* correction for z-axis */
            cpxyz(z, orien[i] + ioffset + 6);
        }

    /*  y-, z- and origin for each base-pair */
    for (i = 1; i <= num_bp; i++) {  /* c6_c8 atoms exist for sure from above */
        avexyz(xyz[c6_c8[1][i]], xyz[c6_c8[2][i]], bp_org + (i - 1) * 3);
        ddxyz(xyz[c6_c8[2][i]], xyz[c6_c8[1][i]], y);
        vec_norm(y);
        ioffset = (i - 1) * 9;
        cpxyz(y, bp_orien + ioffset + 3);

        rnum = pair_num[1][i];
        ib = seidx[rnum][1];
        ie = seidx[rnum][2];
        cehs_base_atoms(AtomName, ib, ie, &num_batom, batom);
        for (j = 1; j <= num_batom; j++)
            cpxyz(xyz[batom[j]], bxyz[j]);

        rnum2 = pair_num[2][i];
        ib2 = seidx[rnum2][1];
        ie2 = seidx[rnum2][2];
        cehs_base_atoms(AtomName, ib2, ie2, &num_batom2, batom2);
        for (j = num_batom + 1; j <= num_batom + num_batom2; j++)
            cpxyz(xyz[batom2[j - num_batom]], bxyz[j]);

        ls_plane(bxyz, num_batom + num_batom2, z, ppos, &odist, adist);

        /* direction of bp follows base I */
        if (dot(z, &orien[1][ioffset + 6]) < 0)
            negate_xyz(z);

        vec_orth(z, y);  /* correction for z-axis */
        cpxyz(z, bp_orien + ioffset + 6);
    }

    ik = 0;  /* check if to reverse z-axis for Z-DNA */

    if (num_bp == 1)
        init_dvector(dorg, 1, 3, 0.0);

    for (i = 1; i <= num_bp; i++) {
        ioffset = (i - 1) * 3;
        if (i < num_bp)
            ddxyz(bp_org + ioffset, bp_org + ioffset + 3, dorg);
        if (dot(dorg, &bp_orien[(i - 1) * 9 + 6]) < 0.0)
            ik++;
    }

    for (i = 1; i <= num_bp; i++) {
        ioffset = (i - 1) * 9;

        if (ik == num_bp)  /* reverse bp z-axis, Z-DNA */
            negate_xyz(bp_orien + ioffset + 6);
        cross(&bp_orien[ioffset + 3], &bp_orien[ioffset + 6], &bp_orien[ioffset]);  /* get x-axis */

        for (j = 1; j <= ds; j++) {  /* each base */
            if (ik == num_bp)  /* reverse z-axis, Z-DNA */
                negate_xyz(orien[j] + ioffset + 6);
            cross(&orien[j][ioffset + 3], &orien[j][ioffset + 6], &orien[j][ioffset]);  /* get x-axis */
        }
    }

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mfi = dmatrix(1, 3, 1, 3);
    bp_par = dmatrix(1, num_bp, 1, 6);
    for (i = 1; i <= num_bp; i++) {
        refs_right_left(i, orien, org, r1, o1, r2, o2);
        cehs_bppar(r1, o1, r2, o2, bp_par[i], mfi, mfoi);
    }

    step_par = dmatrix(1, nbpm1, 1, 6);
    step_str = cmatrix(1, num_bp, 0, 5);
    bphlx = lvector(1, num_bp);  /* helix break marker */
    ie = istart;
    for (i = istart; i <= nbpm1; i++) {
        if (istep > 0) {  /* continuous 1 to 3; 2 to 4 etc */
            ib = i;
            ie = ib + istep;
        } else {  /* next segment: 1 to 3, 3 to 5 etc */
            ib = ie;
            ie = ib - istep;
        }
        if (ie > num_bp)
            break;
        num_step++;
        sprintf(step_str[num_step], "%c%c/%c%c", bp_seq[1][ib],
                bp_seq[1][ie], bp_seq[2][ie], bp_seq[2][ib]);
        refs_i_j(ib, ie, bp_orien, bp_org, r1, o1, r2, o2);
        bz_check(r1, o1, r2, o2, bz, &bz_junction, &z_step);
        bpstep_par(r1, o1, r2, o2, step_par[num_step], mfi, mfoi);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "Local base-pair parameters\n");
    print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "Local step parameters\n");
    prt_stepstr(step_str, num_step, bphlx, 0, step_par, fp);  /* continuous */

    free_dmatrix(bxyz, 1, 2 * NUM_BASE_ATOMS, 1, 3);
    free_dmatrix(orien, 1, ds, 1, num_bp * 9);
    free_dmatrix(org, 1, ds, 1, num_bp * 3);
    free_dmatrix(bp_par, 1, num_bp, 1, 6);
    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mfi, 1, 3, 1, 3);
    free_dmatrix(step_par, 1, nbpm1, 1, 6);
    free_cmatrix(step_str, 1, num_bp, 0, 5);
    free_lvector(bphlx, 1, num_bp);
}

/* Calculate Schnaap linear global parameters: to replace orginal SCHNAP */
void schnaap_global(long num_bp, long num, char **bp_seq, long **chi, double **xyz,
                    double *bp_orien, double *bp_org, FILE * fp)
{
    char *bstr = "      --- ", *fmt = "%10.2f";
    double hrise, std_rise, TipInc, phi;
    double haxis[4], hstart[4], hend[4], hinge[4], dd[7];
    double **rise_twist, **hpars, **r, **rh, **temp, **yhel;
    long ds = 2, i, ioffset, j, k, nbpm1;
    long nvec, C1b[MBASES], C1e[MBASES];
    long **idx;

    idx = lmatrix(1, 4 * num_bp, 1, 2);
    get_CNidx(ds, num_bp, chi, idx, &nvec, C1b, C1e);
    get_axis(nvec, idx, num, xyz, ds, C1b, C1e, &std_rise, &hrise, haxis, hstart, hend);

    if (hrise < EMPTY_CRITERION || std_rise > Gvars.misc_pars.std_curved)
        return;

    r = dmatrix(1, 3, 1, 3);
    temp = dmatrix(1, 3, 1, 3);
    rh = dmatrix(1, 3, 1, 3);
    hpars = dmatrix(1, num_bp, 1, 6);
    yhel = dmatrix(1, num_bp, 1, 3);
    for (i = 1; i <= num_bp; i++) {
        ioffset = (i - 1) * 9;
        orien2mst(bp_orien, ioffset, r);
        TipInc = magang(haxis, &bp_orien[ioffset + 6]);
        cross(haxis, &bp_orien[ioffset + 6], hinge);
        arb_rotation(hinge, -TipInc, temp);
        multi_matrix(temp, 3, 3, r, 3, 3, rh);

        for (j = 1; j <= 3; j++)
            yhel[i][j] = rh[j][2];
        ddxyz(hstart, bp_org + (i - 1) * 3, hend);
        phi = deg2rad(vec_ang(hinge, yhel[i], haxis));
        hpars[i][4] = TipInc * sin(phi);
        hpars[i][5] = TipInc * cos(phi);

        for (j = 1; j <= 3; j++) {
            hpars[i][j] = 0.0;
            for (k = 1; k <= 3; k++)
                hpars[i][j] += hend[k] * rh[k][j];
        }
    }

    nbpm1 = num_bp - 1;
    for (i = 1; i < num_bp; i++) {
        j = i + 1;
        hpars[i][3] = hpars[j][3] - hpars[i][3];
        hpars[i][6] = vec_ang(yhel[i], yhel[j], haxis);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "SCHNAaP global helical parameters\n");
    fprintf(fp, "    step       X-disp    Y-disp   h-Rise     Incl.       Tip" "   h-Twist\n");
    for (i = 1; i <= num_bp; i++) {
        fprintf(fp, " %4ld %c%c%c ", i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);
        for (j = 1; j <= 6; j++)
            if (i < num_bp)
                fprintf(fp, fmt, hpars[i][j]);
            else
                fprintf(fp, (j == 3 || j == 6) ? bstr : fmt, hpars[i][j]);
        fprintf(fp, "\n");
    }
    rise_twist = dmatrix(1, nbpm1, 1, 2);
    if (num_bp > 2) {
        for (i = 1; i < num_bp; i++) {
            rise_twist[i][1] = hpars[i][3];
            rise_twist[i][2] = hpars[i][6];
        }
        ave_dmatrix(hpars, num_bp, 6, dd);
        ave_dmatrix(rise_twist, nbpm1, 2, hend);
        fprintf(fp, "          ");
        print_sep(fp, '~', 60);
        fprintf(fp, "      ave.");
        for (i = 1; i <= 6; i++)
            fprintf(fp, fmt, (i == 3) ? hend[1] : (i == 6) ? hend[2] : dd[i]);
        fprintf(fp, "\n");

        std_dmatrix(hpars, num_bp, 6, dd);
        std_dmatrix(rise_twist, nbpm1, 2, hend);
        fprintf(fp, "      s.d.");
        for (i = 1; i <= 6; i++)
            fprintf(fp, fmt, (i == 3) ? hend[1] : (i == 6) ? hend[2] : dd[i]);
        fprintf(fp, "\n");
    }
    free_lmatrix(idx, 1, 4 * num_bp, 1, 2);
    free_dmatrix(r, 1, 3, 1, 3);
    free_dmatrix(temp, 1, 3, 1, 3);
    free_dmatrix(rh, 1, 3, 1, 3);
    free_dmatrix(hpars, 1, num_bp, 1, 6);
    free_dmatrix(yhel, 1, num_bp, 1, 3);
    free_dmatrix(rise_twist, 1, nbpm1, 1, 2);
}

/* In CEHS: for base-pair parameters, propeller is applied first */
void out_cehs(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org, FILE * fp)
{
    double mfoi[4], o1[4], o2[4];
    double *bp_org, *bp_orien, **bp_par, **mfi, **r1, **r2, **step_par;
    long i, nbpm1;

    nbpm1 = num_bp - 1;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mfi = dmatrix(1, 3, 1, 3);

    bp_org = dvector(1, num_bp * 3);
    bp_orien = dvector(1, num_bp * 9);

    bp_par = dmatrix(1, num_bp, 1, 6);
    step_par = dmatrix(1, nbpm1, 1, 6);

    for (i = 1; i <= num_bp; i++) {
        refs_right_left(i, orien, org, r1, o1, r2, o2);
        cehs_bppar(r1, o1, r2, o2, bp_par[i], mfi, mfoi);

        cpxyz(mfoi, bp_org + (i - 1) * 3);
        mst2orien(bp_orien, (i - 1) * 9, mfi);
    }

    for (i = 1; i <= nbpm1; i++) {
        refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
        bpstep_par(r1, o1, r2, o2, step_par[i], mfi, mfoi);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "CEHS base-pair parameters\n");
    print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "CEHS base-pair step parameters\n");
    prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mfi, 1, 3, 1, 3);
    free_dvector(bp_org, 1, num_bp * 3);
    free_dvector(bp_orien, 1, num_bp * 9);
    free_dmatrix(bp_par, 1, num_bp, 1, 6);
    free_dmatrix(step_par, 1, nbpm1, 1, 6);
}

/* Calculate parameters based on Gorin's scheme */
void compdna(double **rot1, double *org1, double **rot2, double *org2, double *pars,
             double **mst_orien, double *mst_org)
{
    double dorg[4], dz[4], y1_osx[4], y2[4], z1[4], z2[4];
    double xm[4], ym[4], zm[4];
    long i;

    /* use only z- and y-axes for constructing middle frame */
    for (i = 1; i <= 3; i++) {
        z1[i] = rot1[i][3];
        y1_osx[i] = rot1[i][2];
        z2[i] = rot2[i][3];
        y2[i] = rot2[i][2];
    }
    avexyz(org1, org2, mst_org);
    ddxyz(org1, org2, dorg);
    sumxyz(z1, z2, zm);
    ddxyz(z1, z2, dz);
    vec_norm(zm);

    vec_orth(y1_osx, zm);  /* orthogonal y-component */
    vec_orth(y2, zm);
    sumxyz(y1_osx, y2, ym);
    vec_norm(ym);

    cross(ym, zm, xm);

    pars[1] = dot(dorg, xm);
    pars[2] = dot(dorg, ym);
    pars[3] = dot(dorg, zm);
    pars[4] = -2 * rad2deg(asin(dot(dz, ym) / 2));
    pars[5] = 2 * rad2deg(asin(dot(dz, xm) / 2));
    pars[6] = vec_ang(y1_osx, y2, zm);

    x_y_z_2_mtx(xm, ym, zm, mst_orien);
}

void out_compdna(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org, FILE * fp)
{
    double mfoi[4], o1[4], o2[4];
    double *bp_org, *bp_orien, **bp_par, **mfi, **r1, **r2, **step_par;
    long i, nbpm1;

    nbpm1 = num_bp - 1;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mfi = dmatrix(1, 3, 1, 3);

    bp_org = dvector(1, num_bp * 3);
    bp_orien = dvector(1, num_bp * 9);

    bp_par = dmatrix(1, num_bp, 1, 6);
    step_par = dmatrix(1, nbpm1, 1, 6);

    for (i = 1; i <= num_bp; i++) {
        refs_right_left(i, orien, org, r1, o1, r2, o2);
        compdna(r1, o1, r2, o2, bp_par[i], mfi, mfoi);

        cpxyz(mfoi, bp_org + (i - 1) * 3);
        mst2orien(bp_orien, (i - 1) * 9, mfi);
    }

    for (i = 1; i <= nbpm1; i++) {
        refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
        compdna(r1, o1, r2, o2, step_par[i], mfi, mfoi);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "CompDNA base-pair parameters\n");
    print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "CompDNA base-pair step parameters\n");
    prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mfi, 1, 3, 1, 3);
    free_dvector(bp_org, 1, num_bp * 3);
    free_dvector(bp_orien, 1, num_bp * 9);
    free_dmatrix(bp_par, 1, num_bp, 1, 6);
    free_dmatrix(step_par, 1, nbpm1, 1, 6);
}

/* Calculate Curves parameters */
void my_curves(double **rot1, double *org1, double **rot2, double *org2, double *pars)
{
    double dl, du, twist_l, twist_u;
    double j1_osx[4], j2[4], k1[4], k2[4], l1[4], l2[4];
    double d[4], f[4], n[4], pl[4], pu[4], q[4];
    double temp[4], rtmp[4];
    long i;

    /* decompose rot1 and rot2 into JKL */
    for (i = 1; i <= 3; i++) {
        j1_osx[i] = rot1[i][1];
        k1[i] = rot1[i][2];
        l1[i] = rot1[i][3];
        j2[i] = rot2[i][1];
        k2[i] = rot2[i][2];
        l2[i] = rot2[i][3];
    }
    avexyz(org1, org2, q);  /* mean origin */
    sumxyz(l1, l2, n);  /* n is z-axis */
    sumxyz(j1_osx, j2, d);  /* d is x-axis */
    vec_norm(n);
    vec_orth(d, n);
    cross(n, d, f);  /* f is y-axis */

    /* get the intersection of l1 & l2 with the above mean-plane */
    ddxyz(org1, q, temp);
    dl = dot(n, temp) / dot(n, l1);  /* org1 to mean-plane along l1 */
    for (i = 1; i <= 3; i++)
        pl[i] = org1[i] + dl * l1[i];  /* l1 intersection with mean-plane */

    ddxyz(q, org2, temp);
    du = dot(n, temp) / dot(n, l2);  /* org2 to mean-plane along l2 */
    for (i = 1; i <= 3; i++)
        pu[i] = org2[i] - du * l2[i];  /* l2 intersection with mean-plane */

    pars[3] = dl + du;  /* this is why RISE is bigger */
    ddxyz(pl, pu, temp);  /* vector pl ---> pu */
    pars[1] = dot(temp, d);  /* Shift */
    pars[2] = dot(temp, f);  /* Slide */

    cross(l2, d, temp);
    pars[4] = 2 * vec_ang(f, temp, d);  /* Tilt */

    cross(d, temp, rtmp);
    pars[5] = 2 * vec_ang(rtmp, l2, temp);  /* Roll */

    get_vector(f, d, -pars[4] / 2, temp);
    twist_l = vec_ang(k1, temp, l1);

    get_vector(f, d, +pars[4] / 2, temp);
    twist_u = vec_ang(temp, k2, l2);

    pars[6] = twist_l + twist_u;  /* Twist */
}

/* Get mean base-pair frame using Curves method */
void curves_mbt(long ibp, double **orien, double **org, double **cvr, double *cvo)
{
    double o1[4], o2[4], xm[4], ym[4], zm[4];
    double **r1, **r2;
    long i;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);

    refs_right_left(ibp, orien, org, r1, o1, r2, o2);
    for (i = 1; i <= 3; i++) {
        zm[i] = r1[i][3] + r2[i][3];
        xm[i] = r1[i][1] + r2[i][1];
    }
    vec_norm(zm);
    vec_norm(xm);
    cross(zm, xm, ym);  /* xm & zm not orthogonal */
    x_y_z_2_mtx(xm, ym, zm, cvr);
    avexyz(o1, o2, cvo);

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
}

void out_curves(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org, FILE * fp)
{
    double o1[4], o2[4];
    double **bp_par, **mfi, **r1, **r2, **step_par;
    long i, nbpm1;

    nbpm1 = num_bp - 1;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mfi = dmatrix(1, 3, 1, 3);

    bp_par = dmatrix(1, num_bp, 1, 6);
    step_par = dmatrix(1, nbpm1, 1, 6);

    for (i = 1; i <= num_bp; i++) {
        refs_right_left(i, orien, org, r1, o1, r2, o2);
        my_curves(r1, o1, r2, o2, bp_par[i]);
    }

    for (i = 1; i <= nbpm1; i++) {
        curves_mbt(i, orien, org, r1, o1);
        curves_mbt(i + 1, orien, org, r2, o2);
        my_curves(r1, o1, r2, o2, step_par[i]);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "Curves base-pair parameters (by Xiang-Jun Lu)\n");
    print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "Curves base-pair step parameters\n");
    prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mfi, 1, 3, 1, 3);
    free_dmatrix(bp_par, 1, num_bp, 1, 6);
    free_dmatrix(step_par, 1, nbpm1, 1, 6);
}

/* Calculate parameters based on Dickerson's scheme */
void freehelix(double **rot1, double *org1, double **rot2, double *org2, double *pars,
               double **mst_orien, double *mst_org)
{
    double dorg[4], y1_osx[4], y2[4], z1[4], z2[4];
    double xm[4], ym[4], zm[4];
    long i;

    /* use only y- and z-axes for constructing middle frame */
    for (i = 1; i <= 3; i++) {
        y1_osx[i] = rot1[i][2];
        z1[i] = rot1[i][3];
        y2[i] = rot2[i][2];
        z2[i] = rot2[i][3];
    }
    avexyz(org1, org2, mst_org);
    ddxyz(org1, org2, dorg);
    sumxyz(y1_osx, y2, ym);
    sumxyz(z1, z2, zm);
    vec_norm(ym);

    vec_norm(zm);  /* mst-z */
    cross(ym, zm, xm);
    vec_norm(xm);  /* mst-x */
    cross(zm, xm, ym);  /* mst-y */
    vec_norm(ym);

    pars[1] = dot(dorg, xm);
    pars[2] = dot(dorg, ym);
    pars[3] = dot(dorg, zm);
    pars[4] = vec_ang(z1, z2, xm);
    pars[5] = vec_ang(z1, z2, ym);
    pars[6] = vec_ang(y1_osx, y2, zm);

    x_y_z_2_mtx(xm, ym, zm, mst_orien);
}

void out_freehelix(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org,
                   FILE * fp)
{
    double mfoi[4], o1[4], o2[4];
    double *bp_org, *bp_orien, **bp_par, **mfi, **r1, **r2, **step_par;
    long i, nbpm1;

    nbpm1 = num_bp - 1;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mfi = dmatrix(1, 3, 1, 3);

    bp_org = dvector(1, num_bp * 3);
    bp_orien = dvector(1, num_bp * 9);

    bp_par = dmatrix(1, num_bp, 1, 6);
    step_par = dmatrix(1, nbpm1, 1, 6);

    for (i = 1; i <= num_bp; i++) {
        refs_right_left(i, orien, org, r1, o1, r2, o2);
        freehelix(r1, o1, r2, o2, bp_par[i], mfi, mfoi);

        cpxyz(mfoi, bp_org + (i - 1) * 3);
        mst2orien(bp_orien, (i - 1) * 9, mfi);
    }

    for (i = 1; i <= nbpm1; i++) {
        refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
        freehelix(r1, o1, r2, o2, step_par[i], mfi, mfoi);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "FreeHelix base-pair parameters (by Xiang-Jun Lu)\n");
    print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "FreeHelix base-pair step parameters\n");
    prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mfi, 1, 3, 1, 3);
    free_dvector(bp_org, 1, num_bp * 3);
    free_dvector(bp_orien, 1, num_bp * 9);
    free_dmatrix(bp_par, 1, num_bp, 1, 6);
    free_dmatrix(step_par, 1, nbpm1, 1, 6);
}

/* Get the single helical rotation angle and its corresponding axis */
void sgl_helix(double **rot1, double **rot2, double *rot_ang, double *rot_hlx)
{
    double dsum = 0.0, tcos;
    double dx[4] = { EMPTY_NUMBER, 1.0, 0.0, 0.0 };
    double dy[4] = { EMPTY_NUMBER, 0.0, 1.0, 0.0 };
    double **R, **temp;
    long i, ichg = 0, j;

    R = dmatrix(1, 3, 1, 3);
    temp = dmatrix(1, 3, 1, 3);

    transpose_matrix(rot1, 3, 3, temp);
    multi_matrix(temp, 3, 3, rot2, 3, 3, R);  /* rotation matrix w.r.t. 1 */

    for (i = 1; i <= 3; i++)
        dsum += R[i][i];  /* trace of R */
    tcos = 0.5 * (dsum - 1.0);  /* positive rotation angle */
    *rot_ang = (tcos >= 1.0) ? 0.0 : rad2deg(acos(tcos));

    /* helical rotation axis, also = cross(x1-x2, y1-y2) = cross(x2-x1, y2-y1) */
    for (i = 1; i <= 3; i++) {
        dx[i] = R[i][1] - dx[i];
        dy[i] = R[i][2] - dy[i];
    }
    cross(dx, dy, rot_hlx);
    vec_norm(rot_hlx);

    /* check back */
    arb_rotation(rot_hlx, *rot_ang, temp);
    for (i = 1; i <= 3 && !ichg; i++) {
        for (j = 1; j <= 3 && !ichg; j++)
            if (fabs(R[i][j] - temp[i][j]) > XEPS) {
                *rot_ang = -*rot_ang;  /* negate rotation angle */
                ichg = 1;
            }
    }

#if 0  /* September 22, 2006: for verification */
    {
        if (ichg) {
            for (i = 1; i <= 3; i++) {
                for (j = 1; j <= 3; j++)
                    fprintf(stderr, "\t%8.5f(%8.5f)", R[i][j], temp[i][j]);
                fprintf(stderr, "\n");
            }
            print_sep(stderr, '-', 78);
            arb_rotation(rot_hlx, *rot_ang, temp);
            for (i = 1; i <= 3; i++) {
                for (j = 1; j <= 3; j++)
                    fprintf(stderr, "\t%8.5f(%8.5f)", R[i][j], temp[i][j]);
                fprintf(stderr, "\n");
            }
        }
    }
#endif

    free_dmatrix(R, 1, 3, 1, 3);
    free_dmatrix(temp, 1, 3, 1, 3);
}

/* Calculate parameters based on Tung's scheme */
void ngeom(double **rot1, double *org1, double **rot2, double *org2, double *pars,
           double **mst_orien, double *mst_org)
{
    double ang;
    double dorg[4], dorg1[4], haxis[4], tpars[7];
    double **rmtx;
    long i;

    rmtx = dmatrix(1, 3, 1, 3);

    /* translational parameters are the same as in RNA */
    sgl_helix(rot1, rot2, &ang, haxis);

    ddxyz(org1, org2, dorg);
    multi_vec_matrix(dorg, 3, rot1, 3, 3, dorg1);
    arb_rotation(haxis, ang / 2, rmtx);
    multi_vec_matrix(dorg1, 3, rmtx, 3, 3, pars);

    /* rotational parameters are the same as in CEHS */
    bpstep_par(rot1, org1, rot2, org2, tpars, mst_orien, mst_org);

    for (i = 4; i <= 6; i++)
        pars[i] = tpars[i];

    free_dmatrix(rmtx, 1, 3, 1, 3);
}

void out_ngeom(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org, FILE * fp)
{
    double mfoi[4], o1[4], o2[4];
    double *bp_org, *bp_orien, **bp_par, **mfi, **r1, **r2, **step_par;
    long i, nbpm1;

    nbpm1 = num_bp - 1;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mfi = dmatrix(1, 3, 1, 3);

    bp_org = dvector(1, num_bp * 3);
    bp_orien = dvector(1, num_bp * 9);

    bp_par = dmatrix(1, num_bp, 1, 6);
    step_par = dmatrix(1, nbpm1, 1, 6);
    for (i = 1; i <= num_bp; i++) {
        refs_right_left(i, orien, org, r1, o1, r2, o2);
        ngeom(r1, o1, r2, o2, bp_par[i], mfi, mfoi);

        cpxyz(mfoi, bp_org + (i - 1) * 3);
        mst2orien(bp_orien, (i - 1) * 9, mfi);
    }

    for (i = 1; i <= nbpm1; i++) {
        refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
        ngeom(r1, o1, r2, o2, step_par[i], mfi, mfoi);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "NGEOM base-pair parameters\n");
    print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "NGEOM base-pair step parameters\n");
    prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mfi, 1, 3, 1, 3);
    free_dvector(bp_org, 1, num_bp * 3);
    free_dvector(bp_orien, 1, num_bp * 9);
    free_dmatrix(bp_par, 1, num_bp, 1, 6);
    free_dmatrix(step_par, 1, nbpm1, 1, 6);
}

/* Calculate parameters based on Bansal's scheme */
void nuparm(double **rot1, double *org1, double **rot2, double *org2, double *pars,
            double **mst_orien, double *mst_org, double *hpars, long get_hpar)
{
    double a, cx, cy, cz, dx, dy, sx, sy, sz;
    double B[4], dorg[4], x1[4], x2[4], y1_osx[4], y2[4], zh[4];
    double xm[4], ym[4], zm[4];
    double **A, **invA;
    long i;

    for (i = 1; i <= 3; i++) {
        x1[i] = rot1[i][1];
        y1_osx[i] = rot1[i][2];
        x2[i] = rot2[i][1];
        y2[i] = rot2[i][2];
    }
    sumxyz(x1, x2, xm);
    sumxyz(y1_osx, y2, ym);
    avexyz(org1, org2, mst_org);
    ddxyz(org1, org2, dorg);

    /* get the middle frame (xm & ym are not orthogonal) */
    vec_norm(xm);
    vec_norm(ym);
    cross(xm, ym, zm);
    vec_norm(zm);

    x_y_z_2_mtx(xm, ym, zm, mst_orien);

    pars[1] = dot(dorg, xm);
    pars[2] = dot(dorg, ym);
    pars[3] = dot(dorg, zm);
    pars[4] = -2 * rad2deg(asin(dot(y1_osx, zm)));
    pars[5] = 2 * rad2deg(asin(dot(x1, zm)));
    pars[6] = vec_ang(y1_osx, y2, zm);

    if (get_hpar) {  /* helical parameters */
        ddxyz(x2, x1, xm);
        ddxyz(y2, y1_osx, ym);
        cross(xm, ym, zh);
        vec_norm(zh);

        hpars[6] = vec_ang(y1_osx, y2, zh);
        hpars[5] = -rad2deg(asin(dot(zh, x1)));
        hpars[4] = rad2deg(asin(dot(zh, y1_osx)));

        a = deg2rad(hpars[4]);
        cx = cos(a);
        sx = sin(a);
        a = deg2rad(hpars[5]);
        cy = cos(a);
        sy = -sin(a);
        a = deg2rad(pars[6]);  /* not helical twist */
        cz = cos(a);
        sz = sin(a);

        A = dmatrix(1, 3, 1, 3);
        invA = dmatrix(1, 3, 1, 3);

        A[1][1] = 2 * cx * sz;
        A[1][2] = 0.0;
        A[1][3] = 2 * sx;
        A[2][1] = 0.0;
        A[2][2] = -2 * cy * sz;
        A[2][3] = 2 * sy;
        A[3][1] = -2 * cy * sx * sz;
        A[3][2] = 2 * cx * sy * sz;
        A[3][3] = cx * cy * (1 + cz);
        dinverse(A, 3, invA);
        transpose_matrix(invA, 3, 3, A);

        dx = sqrt(2 * (1 + cz + sx * sx * (1 - cz)));
        dy = sqrt(2 * (1 + cz + sy * sy * (1 - cz)));
        B[1] = pars[2] * dy;
        B[2] = pars[1] * dx;
        B[3] = 0.5 * pars[3] * dx * dx;

        multi_vec_matrix(B, 3, A, 3, 3, dorg);
        cpxyz(dorg, hpars);

        free_dmatrix(A, 1, 3, 1, 3);
        free_dmatrix(invA, 1, 3, 1, 3);
    }
}

void out_nuparm(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org, FILE * fp)
{
    double hpi[7], mfoi[4], o1[4], o2[4];
    double *bp_org, *bp_orien, **bp_par, **mfi, **heli_par, **r1, **r2, **step_par;
    long i, nbpm1;

    nbpm1 = num_bp - 1;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mfi = dmatrix(1, 3, 1, 3);

    bp_org = dvector(1, num_bp * 3);
    bp_orien = dvector(1, num_bp * 9);

    bp_par = dmatrix(1, num_bp, 1, 6);
    step_par = dmatrix(1, nbpm1, 1, 6);
    heli_par = dmatrix(1, nbpm1, 1, 6);

    for (i = 1; i <= num_bp; i++) {
        refs_right_left(i, orien, org, r1, o1, r2, o2);
        nuparm(r1, o1, r2, o2, bp_par[i], mfi, mfoi, hpi, 0);

        cpxyz(mfoi, bp_org + (i - 1) * 3);
        mst2orien(bp_orien, (i - 1) * 9, mfi);
    }

    for (i = 1; i <= nbpm1; i++) {
        refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
        nuparm(r1, o1, r2, o2, step_par[i], mfi, mfoi, heli_par[i], 1);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "NUPARM base-pair parameters\n");
    print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "NUPARM base-pair step parameters\n");
    prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "NUPARM base-pair helical parameters\n");
    prt_step_par(bp_seq, num_bp, bphlx, 1, heli_par, fp);

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mfi, 1, 3, 1, 3);
    free_dvector(bp_org, 1, num_bp * 3);
    free_dvector(bp_orien, 1, num_bp * 9);
    free_dmatrix(bp_par, 1, num_bp, 1, 6);
    free_dmatrix(step_par, 1, nbpm1, 1, 6);
    free_dmatrix(heli_par, 1, nbpm1, 1, 6);
}

/* Calculate parameters based on Babcock's scheme */
void rna(double **rot1, double *org1, double *pvt1, double **rot2, double *org2,
         double *pvt2, double *pars, double **mst_orien, double *mst_org)
{
    double ang;
    double dorg[4], dorg1[4], haxis[4], p1[4], p2[4], pt[4];
    double **rmtx;
    long i;

    rmtx = dmatrix(1, 3, 1, 3);

    /* total rotation angle and helical axis */
    sgl_helix(rot1, rot2, &ang, haxis);

    /* rotational parameters */
    for (i = 1; i <= 3; i++)
        pars[i + 3] = ang * haxis[i];
    ddxyz(org1, org2, dorg);

    /* translational parameters */
    multi_vec_matrix(dorg, 3, rot1, 3, 3, dorg1);
    arb_rotation(haxis, ang / 2, rmtx);
    multi_matrix(rot1, 3, 3, rmtx, 3, 3, mst_orien);
    multi_vec_Tmatrix(pvt1, 3, rot1, 3, 3, p1);
    multi_vec_Tmatrix(pvt2, 3, rot2, 3, 3, p2);
    for (i = 1; i <= 3; i++)
        mst_org[i] = 0.5 * (org1[i] + org2[i] + p1[i] + p2[i]);
    ddxyz(pvt1, dorg1, pt);
    multi_vec_matrix(pt, 3, rmtx, 3, 3, p1);
    multi_vec_Tmatrix(pvt2, 3, rmtx, 3, 3, p2);
    for (i = 1; i <= 3; i++)
        pars[i] = p1[i] + p2[i] + pvt1[i] - pvt2[i];

    free_dmatrix(rmtx, 1, 3, 1, 3);
}

/* Make correction to dx and dy using Babcock's pivot point */
void pvt_dxdy(double **rot1, double *org1, double *pvt1, double *pars,
              double **mst_orien, double *mst_org)
{
    double TipInc1, half_rise;
    double axis_h[4], hinge1[4], org1_h[4], t1[4], t2[4];
    double **rot1_h, **temp;
    long i, j;

    temp = dmatrix(1, 3, 1, 3);
    rot1_h = dmatrix(1, 3, 1, 3);

    half_rise = 0.5 * pars[3];
    for (i = 1; i <= 3; i++) {
        t1[i] = rot1[i][3];  /* z1 */
        axis_h[i] = mst_orien[i][3];  /* helical axis */
        org1_h[i] = mst_org[i] - half_rise * mst_orien[i][3];
    }
    TipInc1 = magang(axis_h, t1);
    cross(axis_h, t1, hinge1);
    arb_rotation(hinge1, -TipInc1, temp);
    multi_matrix(temp, 3, 3, rot1, 3, 3, rot1_h);

    for (i = 1; i <= 3; i++) {
        t1[i] = dot(pvt1, rot1[i]);
        negate_xyz(temp[i]);
        temp[i][i] += 1.0;
    }
    for (i = 1; i <= 3; i++)
        t2[i] = dot(t1, temp[i]) + org1[i] - org1_h[i];

    for (i = 1; i <= 2; i++) {
        pars[i] = 0.0;
        for (j = 1; j <= 3; j++)
            pars[i] += t2[j] * rot1_h[j][i];
    }

    free_dmatrix(rot1_h, 1, 3, 1, 3);
    free_dmatrix(temp, 1, 3, 1, 3);
}

void out_rna(long ds, long num_bp, char **bp_seq, long *bphlx, double **orien,
             double **org, FILE * fp)
{
    double pvt0[4] = { EMPTY_NUMBER, 0.0, 0.0, 0.0 };
    double pvt1[4] = { EMPTY_NUMBER, 0.0, 1.808, 0.0 };
    double pvt2[4] = { EMPTY_NUMBER, 0.0, -1.808, 0.0 };
    double hpi[7], mfoi[4], o1[4], o2[4];
    double **mfi, **r1, **r2, **heli_par, **step_par;
    long i, j, k, nbpm1;

    nbpm1 = num_bp - 1;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mfi = dmatrix(1, 3, 1, 3);

    step_par = dmatrix(1, nbpm1, 1, 6);
    heli_par = dmatrix(1, nbpm1, 1, 6);

    if (ds == 2) {  /* double helix */
        double *bp_org, *bp_orien, **bp_par;

        bp_org = dvector(1, num_bp * 3);
        bp_orien = dvector(1, num_bp * 9);
        bp_par = dmatrix(1, num_bp, 1, 6);

        for (i = 1; i <= num_bp; i++) {
            refs_right_left(i, orien, org, r1, o1, r2, o2);
            rna(r1, o1, pvt2, r2, o2, pvt1, bp_par[i], mfi, mfoi);

            cpxyz(mfoi, bp_org + (i - 1) * 3);
            mst2orien(bp_orien, (i - 1) * 9, mfi);
        }

        for (i = 1; i <= nbpm1; i++) {
            refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
            rna(r1, o1, pvt0, r2, o2, pvt0, step_par[i], mfi, mfoi);
            helical_par(r1, o1, r2, o2, heli_par[i], mfi, mfoi);
        }

        print_sep(fp, '*', 76);
        fprintf(fp, "RNA base-pair parameters\n");
        print_par(bp_seq, num_bp, 1, 0, bp_par, fp);

        print_sep(fp, '*', 76);
        fprintf(fp, "RNA base-pair step parameters\n");
        prt_step_par(bp_seq, num_bp, bphlx, 0, step_par, fp);

        print_sep(fp, '*', 76);
        fprintf(fp, "RNA base-pair helical parameters\n");
        prt_step_par(bp_seq, num_bp, bphlx, 1, heli_par, fp);

        free_dvector(bp_org, 1, num_bp * 3);
        free_dvector(bp_orien, 1, num_bp * 9);
        free_dmatrix(bp_par, 1, num_bp, 1, 6);
    }

    /* Step and helical parameters along each strand */
    for (i = 1; i <= ds; i++) {
        for (j = 1; j <= nbpm1; j++) {
            refs_i_j(j, j + 1, orien[i], org[i], r1, o1, r2, o2);
            helical_par(r1, o1, r2, o2, hpi, mfi, mfoi);
            if (i == 1) {
                pvt_dxdy(r1, o1, pvt1, hpi, mfi, mfoi);
                rna(r1, o1, pvt1, r2, o2, pvt1, step_par[j], mfi, mfoi);
            } else {
                pvt_dxdy(r1, o1, pvt2, hpi, mfi, mfoi);
                rna(r1, o1, pvt2, r2, o2, pvt2, step_par[j], mfi, mfoi);
            }
            for (k = 1; k <= 6; k++)
                heli_par[j][k] = hpi[k];
        }
        if (i == 1) {
            print_sep(fp, '*', 76);
            fprintf(fp, "Strand I base step parameters\n");
            print_par(bp_seq, num_bp, 3, 0, step_par, fp);
            print_sep(fp, '*', 76);
            fprintf(fp, "Strand I base helical parameters\n");
            print_par(bp_seq, num_bp, 3, 1, heli_par, fp);
        } else {
            print_sep(fp, '*', 76);
            fprintf(fp, "Strand II base step parameters\n");
            print_par(bp_seq, num_bp, 4, 0, step_par, fp);
            print_sep(fp, '*', 76);
            fprintf(fp, "Strand II base helical parameters\n");
            print_par(bp_seq, num_bp, 4, 1, heli_par, fp);
        }
    }

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mfi, 1, 3, 1, 3);
    free_dmatrix(step_par, 1, nbpm1, 1, 6);
    free_dmatrix(heli_par, 1, nbpm1, 1, 6);
}

/* Calculate parameters for duplex based on all 7 methods */
void other_pars(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org)
{
    FILE *fp;

    fp = open_file(SEVEN_FILE, "w");

    out_cehs(num_bp, bp_seq, bphlx, orien, org, fp);
    out_compdna(num_bp, bp_seq, bphlx, orien, org, fp);
    out_curves(num_bp, bp_seq, bphlx, orien, org, fp);
    out_freehelix(num_bp, bp_seq, bphlx, orien, org, fp);
    out_ngeom(num_bp, bp_seq, bphlx, orien, org, fp);
    out_nuparm(num_bp, bp_seq, bphlx, orien, org, fp);
    out_rna(2L, num_bp, bp_seq, bphlx, orien, org, fp);

    close_file(fp);
}
