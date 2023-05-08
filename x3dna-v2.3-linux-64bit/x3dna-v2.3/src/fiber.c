#include "x3dna.h"

typedef struct {
    long str_id;
    long xml;
    long single;
    long rna;
    long pauling;
    long connect;
    char pdbfile[BUF512];
    char sequence[BUFBIG];
} struct_args;

static void fiber_dir(char *FDIR)
{
    sprintf(FDIR, "%sfiber/", Gvars.X3DNA_HOMEDIR);
}

static void fiber_usage(void)
{
    help3dna_usage("fiber");
}

/* brief description of the 55 fiber structures + Pauling triplex model */
static void list_str_info(char *FDIR)
{
    char README_FILE[BUF512], str[BUF512];
    FILE *fp;

    sprintf(README_FILE, "%sREADME", FDIR);
    fp = open_file(README_FILE, "r");
    while (fgets(str, sizeof str, fp) != NULL)
        if (!strncmp(str, "id# ", 4)) {
            fprintf(stderr, "\n%s", str);
            while (fgets(str, sizeof str, fp) != NULL)
                fprintf(stderr, "%s", str);
        }
    close_file(fp);
    contact_msg(0);
}

static void set_defaults(struct_args * args)
{
    args->str_id = 4;  /* defaulted to B-DNA */
    args->xml = FALSE;
    args->single = FALSE;
    args->rna = FALSE;
    args->pauling = FALSE;
    args->connect = FALSE;
    strcpy(args->pdbfile, "");
    strcpy(args->sequence, "");
}

static void fiber_cmdline(int argc, char *argv[], char *FDIR, struct_args * args)
{
    long i, j;

    if (argc < 2)
        fiber_usage();

    set_defaults(args);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        if (lux_ncmatch(argv[i], "^--?co")) {
            args->connect = TRUE;
            continue;
        }

        if (lux_ncmatch(argv[i], "^--?si")) {
            args->single = TRUE;
            continue;
        }

        if (lux_ncmatch(argv[i], "^--?rn")) {
            args->rna = TRUE;
            continue;
        }

        if (lux_ncmatch(argv[i], "^--?pa")) {
            args->pauling = TRUE;
            if (lux_ncmatch(argv[i], "dna"))
                args->pauling++;
            continue;
        }

        if (sscanf(argv[i], "-%ld", &args->str_id) == 1) {
            if (!lval_in_range(args->str_id, 1, 56)) {
                fprintf(stderr, "Structure id# out of range (1-56)\n");
                fiber_usage();
            } else {
                if (args->str_id == 56)
                    args->pauling = TRUE;
                continue;
            }
        }
        if (lux_ncmatch(argv[i], "^--?se")) {
            get_strvalue(argv[i], args->sequence, FALSE);
            upperstr(args->sequence);
            continue;
        }

        upperstr(argv[i]);

        args->xml = get_xmlArgNumber(argv[i], "-X");
        if (args->xml)
            continue;

        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'M' || argv[i][j] == 'L')
                list_str_info(FDIR);
            else if (argv[i][j] == 'A')
                args->str_id = 1;
            else if (argv[i][j] == 'B')
                args->str_id = 4;
            else if (argv[i][j] == 'C')
                args->str_id = 47;
            else if (argv[i][j] == 'D')
                args->str_id = 48;
            else if (argv[i][j] == 'Z')
                args->str_id = 15;
            else if (argv[i][j] == 'S')
                args->single = TRUE;
            else if (argv[i][j] == 'R')
                args->rna = TRUE;
            else
                fiber_usage();
    }

    if (args->rna)
        args->str_id = 901;  /* based on fiber A-DNA model #1 */

    if (args->pauling)
        args->str_id = 56;  /* Pauling's (wrong) DNA triplex model */

    if (argc == i + 1) {
        char str[BUF512];
        strcpy(args->pdbfile, argv[i]);
        /* sprintf(FDIR, "%sstr%2.2ld/", FDIR, args->str_id); -- problem in fedora 14 64bit version */
        sprintf(str, "%sstr%2.2ld/", FDIR, args->str_id);
        strcpy(FDIR, str);  /* now with fiber structure id included */
    } else
        fiber_usage();

    if (args->single && args->xml) {
        fprintf(stderr, "Currently, the -s option only applies to PDB format\n");
        args->xml = FALSE;
    }
}

static void get_twist_rise(long str_id, double *twist, double *rise)
{
    long idx;
    double twist_rise[55][3] = {
        {1, 32.7, 2.548},
        {2, 65.5, 5.095},
        {3, 0.0, 28.030},
        {4, 36.0, 3.375},
        {5, 72.0, 6.720},
        {6, 180.0, 16.864},
        {7, 38.6, 3.310},
        {8, 40.0, 3.312},
        {9, 120.0, 9.937},
        {10, 80.0, 6.467},
        {11, 80.0, 6.467},
        {12, 45.0, 3.013},
        {13, 90.0, 6.125},
        {14, -90.0, 18.500},
        {15, -60.0, 7.250},
        {16, -51.4, 7.571},
        {17, 0.0, 10.200},
        {18, 36.0, 3.230},
        {19, 36.0, 3.233},
        {20, 32.7, 2.812},
        {21, 30.0, 3.000},
        {22, 32.7, 2.560},
        {23, 32.0, 2.780},
        {24, 36.0, 3.130},
        {25, 32.7, 3.060},
        {26, 36.0, 3.010},
        {27, 32.7, 2.518},
        {28, 32.7, 2.596},
        {29, 32.7, 2.596},
        {30, 32.7, 3.160},
        {31, 30.0, 3.260},
        {32, 32.7, 3.040},
        {33, 30.0, 3.040},
        {34, 30.0, 3.290},
        {35, 31.3, 3.410},
        {36, 60.0, 3.155},
        {37, 36.0, 3.200},
        {38, 36.0, 3.240},
        {39, 72.0, 6.480},
        {40, 72.0, 6.460},
        {41, 144.0, 13.540},
        {42, 32.7, 3.040},
        {43, 36.0, 3.200},
        {44, 36.0, 3.233},
        {45, 36.0, 3.233},
        {46, 36.0, 3.380},
        {47, 40.0, 3.320},
        {48, 87.8, 6.020},
        {49, 60.0, 7.200},
        {50, 60.0, 7.200},
        {51, 31.6, 3.220},
        {52, 90.0, 6.060},
        {53, -38.7, 3.290},
        {54, 32.73, 2.560},
        {55, 36.0, 3.39}
    };

    idx = str_id - 1;
    *twist = twist_rise[idx][1];
    *rise = twist_rise[idx][2];

    fprintf(stderr, "Structure #%ld; ", str_id);
    fprintf(stderr, "Twist: %.1f (degrees); Rise: %.3f (Angstrom)\n", *twist, *rise);
    fflush(stderr);
}

static char *extract_valid_sequence(char *sequence, char *valid_bases, long *num_bp)
{
    long i, k, nb = 0;
    char *bseq;

    k = strlen(sequence);
    bseq = cvector(0, k);

    for (i = 0; i < k; i++) {
        if (is_valid_base(sequence[i], valid_bases))
            bseq[nb++] = sequence[i];
    }
    bseq[nb] = '\0';
    *num_bp = nb;

    return bseq;
}

static void write_coordinates(struct_args * args, long tnum, char **tAtomName,
                              char **tResName, char *tChainID, long *tResSeq, double **txyz)
{
    if (args->xml)
        write_pdbml(args->xml, tnum, tAtomName, tResName, tChainID, tResSeq, txyz, args->pdbfile);
    else {
        if (args->connect)
            atom_lkg(tnum, tAtomName, tResName, tChainID, tResSeq, txyz, FALSE, args->pdbfile);
        else
            write_pdb(tnum, tAtomName, tResName, tChainID, tResSeq, txyz, NULL, args->pdbfile);
    }
}

/* repeating unit is each residue: str# 1, 4, 7, 8 & 12 + (44, 45) + (46, 47) + (53, 55) */
static void residue(double twist, double rise, char *b_rep, long s2r, char *FDIR,
                    struct_args * args)
{
    char fname[BUF512], *Wb, *Cb;
    static char *Wbase = "ACGT", *Cbase = "TGCA";
    static char *Wrbase = "ACGU", *Crbase = "UGCA";
    char *sChainID, *tChainID;
    char **sAtomName, **sResName;
    char **tAtomName, **tAtomName2, **tResName;
    char *bseq, **bp_seq;

    double dphi, dz, phi;
    double **sxyz, **txyz, **txyz2;

    long num_atoms = 0, num_bp;
    long i, ik, j, num1, num2, tnum = 0, tnum2 = 0;
    long idx[4];
    long *sResSeq, *tResSeq;
    long **s2idx;  /* strand II index */

    if (args->rna) {
        Wb = Wrbase;
        Cb = Crbase;
    } else {
        Wb = Wbase;
        Cb = Cbase;
    }

    /* get the base-pair sequence */
    if (is_empty_string(b_rep)) {  /* 1, 4 & 7 + (46, 47) + (53, 55) */
        if (is_empty_string(args->sequence)) {
            bseq = get_sequence(Wb, &num_bp);  /* interactive input */
        } else {  /* specified with -seq command line option */
            bseq = extract_valid_sequence(args->sequence, Wb, &num_bp);
        }
    } else  /* 8 & 12 + (44, 45) */
        bseq = read_repeat(b_rep, 1, Wb, &num_bp);
    bp_seq = single2double(num_bp, bseq, Wb, Cb);  /* bseq is freed */

    /* get the total number of atoms in the final structure */
    for (i = 0; i < 4; i++)
        idx[i] = 0;
    for (i = 1; i <= num_bp; i++)
        for (j = 1; j <= 2; j++)
            ++idx[strchr(Wb, bp_seq[i][j]) - Wb];
    for (i = 0; i < 4; i++) {
        if (idx[i]) {
            sprintf(fname, "%s%c.rpt", FDIR, Wb[i]);
            num_atoms += idx[i] * number_of_atoms(fname, 1, "*");
        }
    }

    tAtomName = cmatrix(1, num_atoms, 0, 4);
    tResName = cmatrix(1, num_atoms, 0, 3);
    tChainID = cvector(1, num_atoms);
    tResSeq = lvector(1, num_atoms);
    txyz = dmatrix(1, num_atoms, 1, 3);
    tAtomName2 = cmatrix(1, num_atoms, 0, 4);
    txyz2 = dmatrix(1, num_atoms, 1, 3);

    sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
    sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
    sChainID = cvector(1, NUM_RESIDUE_ATOMS);
    sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
    sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);

    s2idx = lmatrix(1, num_bp, 1, 2);

    for (i = 1; i <= num_bp; i++) {
        dphi = -(i - 1) * twist;  /* negative */
        dz = -(i - 1) * rise;  /* negative */

        if (s2r == -1) {  /* (46, 47): positive z-axis */
            dphi = -dphi;
            dz = -dz;
        }

        /* Strand I */
        sprintf(fname, "%s%c.rpt", FDIR, bp_seq[i][1]);
        num1 = read_pdb(fname, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL, 1, "*");
        for (j = 1; j <= num1; j++) {
            ik = tnum + j;
            strcpy(tAtomName[ik], sAtomName[j]);
            sprintf(tResName[ik], "  %c", bp_seq[i][1]);
            tChainID[ik] = 'A';
            tResSeq[ik] = i;

            phi = deg2rad(dphi + sxyz[j][3]);
            txyz[ik][1] = sxyz[j][2] * cos(phi);
            txyz[ik][2] = sxyz[j][2] * sin(phi);
            txyz[ik][3] = dz + sxyz[j][1];
        }

        /* Strand II */
        sprintf(fname, "%s%c.rpt", FDIR, bp_seq[i][2]);
        num2 = read_pdb(fname, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL, 1, "*");

        for (j = 1; j <= num2; j++) {
            if (s2r) {  /* 1, 4, 7, 8 & 12 */
                sxyz[j][1] = -sxyz[j][1];  /* -z */
                sxyz[j][3] = -sxyz[j][3];  /* -phi */
            }

            ik = tnum2 + j;
            strcpy(tAtomName2[ik], sAtomName[j]);

            phi = deg2rad(dphi + sxyz[j][3]);
            txyz2[ik][1] = sxyz[j][2] * cos(phi);
            txyz2[ik][2] = sxyz[j][2] * sin(phi);
            txyz2[ik][3] = dz + sxyz[j][1];
        }

        s2idx[i][1] = tnum2;
        s2idx[i][2] = num2;

        tnum += num1;
        tnum2 += num2;
    }

    /* "reverse" strand II and combined with I */
    reverse_stnd2(num_bp, bp_seq, s2idx, &tnum, tAtomName, tResName,
                  tChainID, tResSeq, txyz, tAtomName2, txyz2, 0);

    write_coordinates(args, tnum, tAtomName, tResName, tChainID, tResSeq, txyz);

    free_cmatrix(bp_seq, 1, num_bp, 1, 2);
    free_pdb(num_atoms, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
    free_cmatrix(tAtomName2, 1, num_atoms, 0, 4);
    free_dmatrix(txyz2, 1, num_atoms, 1, 3);
    free_pdb(NUM_RESIDUE_ATOMS, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL);
    free_lmatrix(s2idx, 1, num_bp, 1, 2);
}

static void multi_bp(double twist, double rise, char *FDIR, struct_args * args)
{
    char base_rpt[BUF512], fname[BUF512];
    char *ChainID1, *ChainID2, *tChainID;
    char **AtomName1, **AtomName2, **ResName1, **ResName2, **Miscs1;
    char **tAtomName, **tAtomName2, **tResName, **tResName2;

    double dphi, dz, phi;
    double **xyz1, **xyz2, **txyz, **txyz2;

    long i, ik, j, n_rep, nr, num1, num2, tnum, tnum2;
    long *ResSeq1, *ResSeq2, *tResSeq, *tResSeq2;
    long **s2idx;  /* strand II index */

    switch (args->str_id) {
    case 2:
        strcpy(base_rpt, "ABr5U:ABr5U");
        break;
    case 3:
        strcpy(base_rpt, "ATCGGAATGGT:ACCATTCCGAT");
        break;
    case 5:
        strcpy(base_rpt, "CG:CG");
        break;
    case 6:
        strcpy(base_rpt, "CCCCC:GGGGG");
        break;
    case 9:
        strcpy(base_rpt, "GGT:ACC");
        break;
    case 10:
        strcpy(base_rpt, "AG:CT");
        break;
    case 11:
        strcpy(base_rpt, "AG:CT");
        break;
    case 13:
        strcpy(base_rpt, "CI:CI");
        break;
    case 14:
        strcpy(base_rpt, "ATATAT:ATATAT");
        break;
    case 15:
        strcpy(base_rpt, "GC:GC");
        break;
    case 16:
        strcpy(base_rpt, "As4T:As4T");
        break;
    case 17:
        strcpy(base_rpt, "GC:GC");
        break;
    case 18:
        strcpy(base_rpt, "A:T");
        break;
    case 19:
        strcpy(base_rpt, "A:T");
        break;
    case 20:
        strcpy(base_rpt, "A:U");
        break;
    case 21:
        strcpy(base_rpt, "I:C");
        break;
    case 22:
        strcpy(base_rpt, "A:T");
        break;
    case 23:
        strcpy(base_rpt, "G:C");
        break;
    case 24:
        strcpy(base_rpt, "I:C");
        break;
    case 25:
        strcpy(base_rpt, "A:U");
        break;
    case 26:
        strcpy(base_rpt, "X:X");
        break;
    case 27:
        strcpy(base_rpt, "X:X");
        break;
    case 28:
        strcpy(base_rpt, "s2U:s2U");
        break;
    case 29:
        strcpy(base_rpt, "s2U:s2U");
        break;
    case 37:
        strcpy(base_rpt, "A:U");
        break;
    case 38:
        strcpy(base_rpt, "A:T");
        break;
    case 39:
        strcpy(base_rpt, "AI:CT");
        break;
    case 40:
        strcpy(base_rpt, "AI:CT");
        break;
    case 41:
        strcpy(base_rpt, "AATT:AATT");
        break;
    case 43:
        strcpy(base_rpt, "A:U");
        break;
    case 48:
        strcpy(base_rpt, "AT:AT");
        break;
    case 49:
        strcpy(base_rpt, "CG:CG");
        break;
    case 50:
        strcpy(base_rpt, "GC:GC");
        break;
    case 51:
        strcpy(base_rpt, "A:T");
        break;
    case 52:
        strcpy(base_rpt, "AT:AT");
        break;
    default:
        fatal("wrong structure id# for this function\n");
    }
    fprintf(stderr, "Repeating unit: %s\n", base_rpt);

    n_rep = repeat_num();

    sprintf(fname, "%sS1.rpt", FDIR);
    num1 = number_of_atoms(fname, 1, "*");
    AtomName1 = cmatrix(1, num1, 0, 4);
    ResName1 = cmatrix(1, num1, 0, 3);
    ChainID1 = cvector(1, num1);
    ResSeq1 = lvector(1, num1);
    xyz1 = dmatrix(1, num1, 1, 3);
    Miscs1 = cmatrix(1, num1, 0, NMISC);
    read_pdb(fname, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, Miscs1, 1, "*");
    residue_idx(num1, ResSeq1, Miscs1, ChainID1, ResName1, &nr);

    sprintf(fname, "%sS2.rpt", FDIR);
    num2 = number_of_atoms(fname, 1, "*");
    AtomName2 = cmatrix(1, num2, 0, 4);
    ResName2 = cmatrix(1, num2, 0, 3);
    ChainID2 = cvector(1, num2);
    ResSeq2 = lvector(1, num2);
    xyz2 = dmatrix(1, num2, 1, 3);
    read_pdb(fname, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz2, NULL, 1, "*");
    tnum2 = n_rep * num2;

    tnum = n_rep * num1 + tnum2;

    tChainID = cvector(1, tnum);
    tAtomName = cmatrix(1, tnum, 0, 4);
    tResName = cmatrix(1, tnum, 0, 3);
    tResSeq = lvector(1, tnum);
    txyz = dmatrix(1, tnum, 1, 3);

    tAtomName2 = cmatrix(1, tnum2, 0, 4);
    tResName2 = cmatrix(1, tnum2, 0, 3);
    tResSeq2 = lvector(1, tnum2);
    txyz2 = dmatrix(1, tnum2, 1, 3);

    s2idx = lmatrix(1, n_rep, 1, 2);

    tnum = 0;
    tnum2 = 0;
    for (i = 1; i <= n_rep; i++) {
        if (lval_in_range(args->str_id, 38, 41) || lval_in_range(args->str_id, 48, 52)) {  /* 38-41/48-52 */
            dphi = (i - 1) * twist;
            dz = (i - 1) * rise;
        } else {
            dphi = -(i - 1) * twist;  /* negative */
            dz = -(i - 1) * rise;  /* negative */
        }

        /* Strand I */
        for (j = 1; j <= num1; j++) {
            ik = tnum + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], ResName1[j]);
            tChainID[ik] = 'A';
            tResSeq[ik] = ResSeq1[j] + (i - 1) * nr;

            phi = deg2rad(dphi + xyz1[j][3]);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        /* Strand II */
        for (j = 1; j <= num2; j++) {
            ik = tnum2 + j;
            strcpy(tAtomName2[ik], AtomName2[j]);
            strcpy(tResName2[ik], ResName2[j]);
            tResSeq2[ik] = ResSeq2[j] + nr * n_rep + (n_rep - i) * nr;

            phi = deg2rad(dphi + xyz2[j][3]);
            txyz2[ik][1] = xyz2[j][2] * cos(phi);
            txyz2[ik][2] = xyz2[j][2] * sin(phi);
            txyz2[ik][3] = dz + xyz2[j][1];
        }

        s2idx[i][1] = tnum2;
        s2idx[i][2] = num2;

        tnum += num1;
        tnum2 += num2;
    }

    /* "reverse" strand II and combined with I */
    for (i = n_rep; i >= 1; i--) {
        for (j = s2idx[i][1] + 1; j <= s2idx[i][1] + s2idx[i][2]; j++) {
            ik = tnum + j - s2idx[i][1];
            strcpy(tAtomName[ik], tAtomName2[j]);
            strcpy(tResName[ik], tResName2[j]);
            tChainID[ik] = 'B';
            tResSeq[ik] = tResSeq2[j];
            cpxyz(txyz2[j], txyz[ik]);
        }
        tnum += s2idx[i][2];
    }

    write_coordinates(args, tnum, tAtomName, tResName, tChainID, tResSeq, txyz);

    free_pdb(num1, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, Miscs1);
    free_pdb(num2, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz2, NULL);
    free_pdb(tnum, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
    free_pdb(tnum2, NULL, tAtomName2, tResName2, NULL, tResSeq2, txyz2, NULL);
    free_lmatrix(s2idx, 1, n_rep, 1, 2);
}

/* structures 30-34 & 42 */
static void triplex(double twist, double rise, char *FDIR, struct_args * args)
{
    char triplet[BUF512], fname[BUF512];
    char *ChainID1, *ChainID2, *tChainID;
    char **AtomName1, **AtomName2, **ResName1, **ResName2;
    char **tAtomName, **tAtomName2, **tResName, **tResName2;
    char **bp_seq;  /* for calling reverse_stnd2 only */

    double dphi, dz, phi;
    double **xyz1, **xyz2, **xyz3, **txyz, **txyz2;

    long i, ik, j, n_rep, num1, num2, tnum, tnum2, offset3;
    long *ResSeq1, *ResSeq2, *tResSeq, *tResSeq2;
    long **s2idx;  /* strand II index */

    switch (args->str_id) {
    case 30:
        strcpy(triplet, "CIC");
        break;
    case 31:
        strcpy(triplet, "TAT");
        break;
    case 32:
        strcpy(triplet, "UAU");
        break;
    case 33:
        strcpy(triplet, "UAU");
        break;
    case 34:
        strcpy(triplet, "IAI");
        break;
    case 42:
        strcpy(triplet, "UAU");
        break;
    default:
        fatal("wrong structure id# for function <triplex>\n");
        break;
    }  /* strands 1 & 3 have the same residue */
    fprintf(stderr, "Triplex repeating unit: %s\n", triplet);

    n_rep = repeat_num();

    /* I: middle strand, R residue, 5'-->3' */
    sprintf(fname, "%s%c2.rpt", FDIR, triplet[1]);
    num1 = number_of_atoms(fname, 1, "*");
    AtomName1 = cmatrix(1, num1, 0, 4);
    ResName1 = cmatrix(1, num1, 0, 3);
    ChainID1 = cvector(1, num1);
    ResSeq1 = lvector(1, num1);
    xyz1 = dmatrix(1, num1, 1, 3);
    read_pdb(fname, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL, 1, "*");

    /* II: 3'-->5', anti-parallel to the middle strand */
    sprintf(fname, "%s%c1.rpt", FDIR, triplet[0]);
    num2 = number_of_atoms(fname, 1, "*");
    AtomName2 = cmatrix(1, num2, 0, 4);
    ResName2 = cmatrix(1, num2, 0, 3);
    ChainID2 = cvector(1, num2);
    ResSeq2 = lvector(1, num2);
    xyz2 = dmatrix(1, num2, 1, 3);
    read_pdb(fname, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz2, NULL, 1, "*");
    tnum2 = n_rep * num2;

    /* III: 5'-->3' parallel to the middle strand */
    sprintf(fname, "%s%c3.rpt", FDIR, triplet[2]);
    xyz3 = dmatrix(1, num2, 1, 3);
    read_pdb(fname, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz3, NULL, 1, "*");

    offset3 = n_rep * num1 + tnum2;
    tnum = offset3 + tnum2;

    tChainID = cvector(1, tnum);
    tAtomName = cmatrix(1, tnum, 0, 4);
    tResName = cmatrix(1, tnum, 0, 3);
    tResSeq = lvector(1, tnum);
    txyz = dmatrix(1, tnum, 1, 3);

    tAtomName2 = cmatrix(1, tnum2, 0, 4);
    tResName2 = cmatrix(1, tnum2, 0, 3);
    tResSeq2 = lvector(1, tnum2);
    txyz2 = dmatrix(1, tnum2, 1, 3);

    s2idx = lmatrix(1, n_rep, 1, 2);

    tnum = 0;
    tnum2 = 0;
    for (i = 1; i <= n_rep; i++) {
        dphi = -(i - 1) * twist;  /* negative */
        dz = -(i - 1) * rise;  /* negative */

        /* Strand I */
        for (j = 1; j <= num1; j++) {
            ik = tnum + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            sprintf(tResName[ik], "  %c", triplet[1]);
            tChainID[ik] = 'A';
            tResSeq[ik] = i;

            phi = deg2rad(dphi + xyz1[j][3]);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        /* Strand II */
        for (j = 1; j <= num2; j++) {
            ik = tnum2 + j;
            strcpy(tAtomName2[ik], AtomName2[j]);

            phi = deg2rad(dphi + xyz2[j][3]);
            txyz2[ik][1] = xyz2[j][2] * cos(phi);
            txyz2[ik][2] = xyz2[j][2] * sin(phi);
            txyz2[ik][3] = dz + xyz2[j][1];
        }

        /* Strand III */
        for (j = 1; j <= num2; j++) {
            ik = offset3 + j;
            strcpy(tAtomName[ik], AtomName2[j]);
            sprintf(tResName[ik], "  %c", triplet[2]);
            tChainID[ik] = 'C';
            tResSeq[ik] = 2 * n_rep + i;

            phi = deg2rad(dphi + xyz3[j][3]);
            txyz[ik][1] = xyz3[j][2] * cos(phi);
            txyz[ik][2] = xyz3[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz3[j][1];
        }

        s2idx[i][1] = tnum2;
        s2idx[i][2] = num2;

        tnum += num1;
        tnum2 += num2;
        offset3 += num2;
    }

    /* "reverse" strand II and combined with I */
    bp_seq = cmatrix(1, n_rep, 1, 2);
    for (i = 1; i <= n_rep; i++)
        bp_seq[i][2] = triplet[0];
    reverse_stnd2(n_rep, bp_seq, s2idx, &tnum, tAtomName, tResName,
                  tChainID, tResSeq, txyz, tAtomName2, txyz2, 0);
    tnum += tnum2;

    write_coordinates(args, tnum, tAtomName, tResName, tChainID, tResSeq, txyz);

    free_pdb(num1, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL);
    free_pdb(num2, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz2, NULL);
    free_dmatrix(xyz3, 1, num2, 1, 3);
    free_pdb(tnum, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
    free_pdb(tnum2, NULL, tAtomName2, tResName2, NULL, tResSeq2, txyz2, NULL);
    free_lmatrix(s2idx, 1, n_rep, 1, 2);
    free_cmatrix(bp_seq, 1, n_rep, 1, 2);
}

static void extract_sequence_to_3chains(char *sequence, char *seq1, char *seq2, char *seq3)
{
    char *seqOK, *items[4], *C10 = "CCCCCCCCCC", bases[BUF32] = "ACGTU";
    long len, nitem, nitem_max = 3;

    len = strlen(sequence);
    if (!len) {  /* just --pauling, without --sequence */
        strcpy(seq1, C10);
        strcpy(seq2, C10);
        strcpy(seq3, C10);
        return;
    }

    /* --sequence=AAATTT; GGUUUU ACGTT --- check up to three */
    nitem = item_list(sequence, items, nitem_max, ",;|: ");
    if (!nitem)
        fatal("[e] invalid --sequence option\n");

    seqOK = extract_valid_sequence(items[1], bases, &len);
    assert(len && len < BUF2K);
    strcpy(seq1, seqOK);

    strcpy(seq2, seqOK);
    strcpy(seq3, seqOK);

    if (nitem == 1) {
        free_cvector(seqOK, 0, DUMMY);
        return;
    }

    free_cvector(seqOK, 0, DUMMY);

    strcat(bases, "0");  /* special case to remove 2nd/3rd chain */

    seqOK = extract_valid_sequence(items[2], bases, &len);
    assert(len < BUF2K);
    strcpy(seq2, seqOK);

    if (nitem == 2) {
        free_cvector(seqOK, 0, DUMMY);
        return;
    }

    free_cvector(seqOK, 0, DUMMY);

    seqOK = extract_valid_sequence(items[3], bases, &len);
    assert(len && len < BUF2K);
    strcpy(seq3, seqOK);

    free_cvector(seqOK, 0, DUMMY);
}

static void check_null_strand(char *seq)
{
    if (strchr(seq, '0'))
        strcpy(seq, "");
}

static long get_numAtoms_perChain(char *FDIR, char *seq)
{
    char fname[BUF512];
    long i, num = 0;

    for (i = 0; i < (long) strlen(seq); i++) {
        sprintf(fname, "%s%c.rpt", FDIR, seq[i]);
        num += number_of_atoms(fname, 1, "*");
    }

    return num;
}

static void populate_pdb_perChain(long chain_number, char *seq, long *idx, long dna,
                                  char *FDIR, char **tAtomName, char **tResName,
                                  char *tChainID, long *tResSeq, double **txyz)
{
    char *chains = ".ABC";
    char base, fname[BUF512], *sChainID, **sAtomName, **sResName;
    double rho, phi, z, d_phi, d_z, **sxyz, chain_phi[] = { -1, 0, 120, 240 };
    long i, j, num, tnum = *idx, *sResSeq;
    long nbases = strlen(seq);

    if (!nbases)
        return;

    sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
    sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
    sChainID = cvector(1, NUM_RESIDUE_ATOMS);
    sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
    sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);

    for (i = 0; i < nbases; i++) {
        base = seq[i];
        sprintf(fname, "%s%c.rpt", FDIR, base);
        num = read_pdb(fname, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL, 1, "*");

        d_z = i * 3.4;
        d_phi = i * 105 + chain_phi[chain_number];

        for (j = 1; j <= num; j++) {
            if (dna && is_equal_string(sAtomName[j], " O2'"))
                continue;

            tnum++;

            strcpy(tAtomName[tnum], sAtomName[j]);
            sprintf(tResName[tnum], "  %c", base);
            tChainID[tnum] = chains[chain_number];
            tResSeq[tnum] = i + 1;

            z = sxyz[j][1] + d_z;
            rho = sxyz[j][2];
            phi = deg2rad(sxyz[j][3] + d_phi);

            txyz[tnum][1] = rho * cos(phi);
            txyz[tnum][2] = rho * sin(phi);
            txyz[tnum][3] = z;
        }
    }

    free_pdb(NUM_RESIDUE_ATOMS, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL);

    *idx = tnum;
}

/* Pauling, L., & Corey, R. B. (1953). A proposed structure for the
 * nucleic acids. Proceedings of the National Academy of Sciences,
 * 39(2), 84-97. --- The model is derived based on Table 1 */
static void pauling_triplex(char *FDIR, struct_args * args)
{
    char seq1[BUF2K], seq2[BUF2K], seq3[BUF2K];
    char *tChainID, **tAtomName, **tResName;
    double **txyz;
    long num_atoms, *tResSeq, k = 0, dna = args->pauling > TRUE;
    long n1, n2, n3;

    extract_sequence_to_3chains(args->sequence, seq1, seq2, seq3);
    check_null_strand(seq2);
    check_null_strand(seq3);

    n1 = strlen(seq1);
    n2 = strlen(seq2);
    n3 = strlen(seq3);
    if (n1 != n2 || n1 != n3)
        fprintf(stderr, "[w] chains differ in lengths: %ld vs %ld vs %ld\n", n1, n2, n3);

    num_atoms = get_numAtoms_perChain(FDIR, seq1) + get_numAtoms_perChain(FDIR, seq2) +
        get_numAtoms_perChain(FDIR, seq3);

    tAtomName = cmatrix(1, num_atoms, 0, 4);
    tResName = cmatrix(1, num_atoms, 0, 3);
    tChainID = cvector(1, num_atoms);
    tResSeq = lvector(1, num_atoms);
    txyz = dmatrix(1, num_atoms, 1, 3);

    populate_pdb_perChain(1, seq1, &k, dna, FDIR, tAtomName, tResName, tChainID, tResSeq, txyz);
    populate_pdb_perChain(2, seq2, &k, dna, FDIR, tAtomName, tResName, tChainID, tResSeq, txyz);
    populate_pdb_perChain(3, seq3, &k, dna, FDIR, tAtomName, tResName, tChainID, tResSeq, txyz);

    write_coordinates(args, k, tAtomName, tResName, tChainID, tResSeq, txyz);

    free_pdb(num_atoms, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
}

/* structure 35 Poly(I).Poly(I).Poly(I).Poly(I) */
static void quadruplex(double twist, double rise, char *FDIR, struct_args * args)
{
    char fname[BUF512];
    char *ChainID1, *tChainID;
    char **AtomName1, **ResName1, **tAtomName, **tResName;
    double dphi, dz, phi;
    double **xyz1, **txyz;
    long i, ik, j, n_rep, num1, tnum;
    long offset1, offset2, offset3, offset4;
    long *ResSeq1, *tResSeq;

    /* all 4 strands are parallel */
    fprintf(stderr, "Quadruplex repeating unit: I.I.I.I\n");

    n_rep = repeat_num();

    sprintf(fname, "%sI.rpt", FDIR);
    num1 = number_of_atoms(fname, 1, "*");
    AtomName1 = cmatrix(1, num1, 0, 4);
    ResName1 = cmatrix(1, num1, 0, 3);
    ChainID1 = cvector(1, num1);
    ResSeq1 = lvector(1, num1);
    xyz1 = dmatrix(1, num1, 1, 3);
    read_pdb(fname, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL, 1, "*");

    offset1 = 0;
    offset2 = n_rep * num1;
    offset3 = 2 * n_rep * num1;
    offset4 = 3 * n_rep * num1;
    tnum = 4 * n_rep * num1;

    tChainID = cvector(1, tnum);
    tAtomName = cmatrix(1, tnum, 0, 4);
    tResName = cmatrix(1, tnum, 0, 3);
    tResSeq = lvector(1, tnum);
    txyz = dmatrix(1, tnum, 1, 3);

    for (i = 1; i <= n_rep; i++) {
        dphi = -(i - 1) * twist;  /* negative */
        dz = -(i - 1) * rise;  /* negative */

        /* Strand I */
        for (j = 1; j <= num1; j++) {
            ik = offset1 + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], "  I");
            tChainID[ik] = 'A';
            tResSeq[ik] = i;

            phi = deg2rad(dphi + xyz1[j][3]);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        /* Strand II */
        for (j = 1; j <= num1; j++) {
            ik = offset2 + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], "  I");
            tChainID[ik] = 'B';
            tResSeq[ik] = n_rep + i;

            phi = deg2rad(dphi + xyz1[j][3] + 90.0);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        /* Strand III */
        for (j = 1; j <= num1; j++) {
            ik = offset3 + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], "  I");
            tChainID[ik] = 'C';
            tResSeq[ik] = 2 * n_rep + i;

            phi = deg2rad(dphi + xyz1[j][3] + 180.0);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        /* Strand IV */
        for (j = 1; j <= num1; j++) {
            ik = offset4 + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], "  I");
            tChainID[ik] = 'D';
            tResSeq[ik] = 3 * n_rep + i;

            phi = deg2rad(dphi + xyz1[j][3] - 90.0);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        offset1 += num1;
        offset2 += num1;
        offset3 += num1;
        offset4 += num1;
    }

    write_coordinates(args, tnum, tAtomName, tResName, tChainID, tResSeq, txyz);

    free_pdb(num1, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL);
    free_pdb(tnum, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
}

/* structure 36 Poly(eC) */
static void shelix(double twist, double rise, char *FDIR, struct_args * args)
{
    char fname[BUF512];
    char *ChainID1, *tChainID;
    char **AtomName1, **ResName1, **tAtomName, **tResName;
    double dphi, dz, phi;
    double **xyz1, **txyz;
    long i, ik, j, n_rep, num1, tnum;
    long *ResSeq1, *tResSeq;

    fprintf(stderr, "Single helix repeating unit: eC\n");

    n_rep = repeat_num();

    sprintf(fname, "%sC.rpt", FDIR);
    num1 = number_of_atoms(fname, 1, "*");
    AtomName1 = cmatrix(1, num1, 0, 4);
    ResName1 = cmatrix(1, num1, 0, 3);
    ChainID1 = cvector(1, num1);
    ResSeq1 = lvector(1, num1);
    xyz1 = dmatrix(1, num1, 1, 3);
    read_pdb(fname, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL, 1, "*");

    tnum = n_rep * num1;

    tChainID = cvector(1, tnum);
    tAtomName = cmatrix(1, tnum, 0, 4);
    tResName = cmatrix(1, tnum, 0, 3);
    tResSeq = lvector(1, tnum);
    txyz = dmatrix(1, tnum, 1, 3);

    tnum = 0;
    for (i = 1; i <= n_rep; i++) {
        dphi = -(i - 1) * twist;  /* negative */
        dz = -(i - 1) * rise;  /* negative */

        for (j = 1; j <= num1; j++) {
            ik = tnum + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], "  C");
            tChainID[ik] = 'A';
            tResSeq[ik] = i;

            phi = deg2rad(dphi + xyz1[j][3]);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        tnum += num1;
    }

    write_coordinates(args, tnum, tAtomName, tResName, tChainID, tResSeq, txyz);

    free_pdb(num1, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL);
    free_pdb(tnum, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
}

static void extract_single_strand(char *pdbfile, long connect)
{
    char *p, tempfile[BUF512] = "x3dna_v2_temp_fiber.pdb";
    long k = 0, max_snum_A;
    FILE *fp0, *fp1;

    remove_file(tempfile);
    rename_file(pdbfile, tempfile);

    fp0 = open_file(tempfile, "r");
    fp1 = open_file(pdbfile, "w");

    fprintf(fp1, "REMARK    %s\n", Gvars.X3DNA_VER);

    while ((p = my_getline(fp0)) != NULL) {
        if (str_pmatch(p, "ATOM  ") && strlen(p) >= 22 && (p[21] == 'A' || p[21] == 'a')) {
            fprintf(fp1, "%s\n", p);
            if (connect)
                sscanf(p + 6, "%ld", &max_snum_A);
        }
        free(p);
    }

    if (connect) {
        rewind(fp0);
        while ((p = my_getline(fp0)) != NULL) {
            if (str_pmatch(p, "CONECT") && strlen(p) >= 11 &&
                sscanf(p + 6, "%ld", &k) == 1 && k <= max_snum_A)
                fprintf(fp1, "%s\n", p);
            free(p);
        }
    }

    fprintf(fp1, "END\n");

    close_file(fp0);
    close_file(fp1);
}

/* utility for generating 55 fiber DNA/RNA structures */
int main(int argc, char *argv[])
{
    struct_args args;

    char FDIR[BUF512];
    long ABC[] = { -1, 1, 4, 7 };
    long S1S2[] = { -1, 2, 3, 5, 6, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20,
        21, 22, 23, 24, 25, 26, 27, 28, 29, 37, 38, 39, 40, 41, 43,
        48, 49, 50, 51, 52
    };
    long TRI[] = { -1, 30, 31, 32, 33, 34, 42 };
    long str_id;
    double rise, twist;

    set_my_globals(argv[0]);
    fiber_dir(FDIR);

    fiber_cmdline(argc, argv, FDIR, &args);

    if (args.pauling) {
        pauling_triplex(FDIR, &args);
        goto TIDY_UP;
    }

    str_id = args.str_id % 900;
    get_twist_rise(str_id, &twist, &rise);

    ABC[0] = sizeof ABC / sizeof ABC[0] - 1;
    S1S2[0] = sizeof S1S2 / sizeof S1S2[0] - 1;
    TRI[0] = sizeof TRI / sizeof TRI[0] - 1;

    /* Five types of fiber structure generating routines
       % o single residue as repeating unit, including 1, 4, 7, 8 & 12,
       %                                      (44,45), (46,47), (53-55)
       % o base-pair(s) as repeating unit, including 2, 3, 5, 6, 9-11,
       %                                         13-29, 37-41 & 43, 48-52
       % o triplex, including 30-34, 42
       % o quadruplex, including 35
       % o single helix, including 36 */
    if (lval_in_set(str_id, 1, ABC[0], ABC))
        residue(twist, rise, "", 1, FDIR, &args);
    else if (str_id == 8)  /* C-DNA  poly d(GGT) : poly d(ACC) */
        residue(twist, rise, "GGT", 1, FDIR, &args);
    else if (str_id == 12)  /* poly d(AAT) : poly d(ATT) */
        residue(twist, rise, "AAT", 1, FDIR, &args);
    else if (lval_in_range(str_id, 44, 45))  /* poly d(A) : poly d(T) */
        residue(twist, rise, "A", 0, FDIR, &args);  /* T in stnd 2: positive z-axis */
    else if (lval_in_range(str_id, 46, 47) ||  /* BI/BII-DNA: Levitt JMB 2000 */
             lval_in_range(str_id, 53, 55))  /* Premilat & Albiser A-, B- and C-DNA */
        residue(twist, rise, "", -1, FDIR, &args);  /* positive z-axis */
    else if (lval_in_set(str_id, 1, S1S2[0], S1S2))
        multi_bp(twist, rise, FDIR, &args);
    else if (lval_in_set(str_id, 1, TRI[0], TRI))
        triplex(twist, rise, FDIR, &args);
    else if (str_id == 35)
        quadruplex(twist, rise, FDIR, &args);
    else if (str_id == 36)
        shelix(twist, rise, FDIR, &args);
    else
        fatal("fiber structure number [%ld] out of 1-55 range\n", str_id);

  TIDY_UP:
    if (args.single)
        extract_single_strand(args.pdbfile, args.connect);

    clear_my_globals();

    return 0;
}
