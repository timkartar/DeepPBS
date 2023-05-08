#include "x3dna.h"

typedef struct {
    char inpfile[BUF512];
    char outfile[BUF512];
    char pdb[BUF512];  /* convert to PDB format */
    long mol;  /* MDL molfile v2000 or v3000 */
    double scale_factor;
} struct_args;

static void cvt2mol(long molVersion, char *inpfile, char *outfile)
{
    char **AtomName;
    double **xyz;
    long i, j, nbond, num, zero = 0, *ibase, **linkage;
    time_t current_time;
    FILE *fp;

    /* read in ALCHEMY file */
    get_alc_nums(inpfile, &num, &nbond);
    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 0, 2);
    read_alc(inpfile, &num, &nbond, AtomName, xyz, ibase, linkage);

    fp = open_file(outfile, "w");
    fprintf(fp, "%s\n", inpfile);
    fprintf(fp, "XL 3DNAv2 \n");
    time(&current_time);
    fprintf(fp, "Converted from Alchemy format: %s", ctime(&current_time));

    if (molVersion == 3 || num > 999 || nbond > 999) {  /* v3000 */
        fprintf(fp, "%3ld%3ld%3ld   %3ld%3ld            999 V3000\n", zero, zero,
                zero, zero, zero);
        fprintf(fp, "M  V30 BEGIN CTAB\n");
        fprintf(fp, "M  V30 COUNTS %5ld %5ld 0 0 0\n", num, nbond);
        fprintf(fp, "M  V30 BEGIN ATOM\n");
        for (i = 1; i <= num; i++)
            fprintf(fp, "M  V30 %5ld %-2s %10.4f %10.4f %10.4f 0\n", i, AtomName[i],
                    xyz[i][1], xyz[i][2], xyz[i][3]);
        fprintf(fp, "M  V30 END ATOM\n");
        fprintf(fp, "M  V30 BEGIN BOND\n");
        for (i = 1; i <= nbond; i++)
            fprintf(fp, "M  V30 %5ld 1 %4ld %4ld\n", i, linkage[i][1], linkage[i][2]);
        fprintf(fp, "M  V30 END BOND\n");
        fprintf(fp, "M  V30 END CTAB\n");
    } else {  /* v2000 */
        fprintf(fp, "%3ld%3ld%3ld   %3ld%3ld              1 V2000\n", num, nbond,
                zero, zero, zero);
        for (i = 1; i <= num; i++) {
            for (j = 1; j <= 3; j++) {
                fprintf(fp, "%10.4f", xyz[i][j]);
            }
            fprintf(fp, " %-3s 0  0  0  0  0  0  0  0  0  0  0  0\n", AtomName[i]);
        }
        for (i = 1; i <= nbond; i++) {
            fprintf(fp, "%3ld%3ld  1  0  0  0\n", linkage[i][1], linkage[i][2]);
        }
    }
    fprintf(fp, "M END\n");

    close_file(fp);

    free_alc(num, nbond, AtomName, xyz, ibase, 0, linkage);
}

static void convert_atom_hetatm(char *pdbinp, char *pdbout)
{
    char *p0;
    FILE *fpi, *fpo;

    fpi = open_file(pdbinp, "r");
    fpo = open_file(pdbout, "w");

    while ((p0 = my_getline(fpi)) != NULL) {
        if (str_pmatch(p0, "ATOM  "))
            strncpy(p0, "HETATM", 6);
        fprintf(fpo, "%s\n", p0);
        free(p0);
    }

    close_file(fpi);
    close_file(fpo);
}

static void cvt2pdb(struct_args * args)
{
    char **AtomAlc, **AtomName, **ResName, *ChainID;
    char rname[BUF32] = "", cid = 'A', *items[BUF32];
    char *alc2pdb = "alc-atom.pdb";
    double **xyzAlc, **xyz;
    long i, nbond, num, nitem, rseq = 1, *ibase, **linkage;
    long *ResSeq, **connect;

    /* processing of -pdb=XXX:Y:ZZZZ option */
    upperstr(args->pdb);
    nitem = item_list(args->pdb, items, 3, "/;:");
    if (nitem > 0)
        strncat(rname, items[1], 3);
    if (nitem > 1) {
        cid = items[2][0];
        if (!isprint((int) cid))
            cid = 'A';
    }
    if (nitem > 2) {
        rseq = cvt2long(items[3]);
        if (rseq == LONG_MAX)
            rseq = 1;
    }

    /* read in ALCHEMY file */
    get_alc_nums(args->inpfile, &num, &nbond);
    AtomAlc = cmatrix(1, num, 0, 2);
    xyzAlc = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 0, 2);
    read_alc(args->inpfile, &num, &nbond, AtomAlc, xyzAlc, ibase, linkage);

    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    for (i = 1; i <= num; i++) {
        sprintf(AtomName[i], " %-2.2s ", AtomAlc[i]);
        sprintf(ResName[i], "%3.3s", rname);
        ChainID[i] = cid;
        ResSeq[i] = rseq;
        cpxyz(xyzAlc[i], xyz[i]);
    }

    connect = lmatrix(1, num, 1, 7);
    lkg2connect(AtomName, 1, num, nbond, linkage, connect);

    write_pdbcnt(num, AtomName, ResName, ChainID, ResSeq, xyz, connect, alc2pdb);
    convert_atom_hetatm(alc2pdb, args->outfile);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, NULL);
    free_lmatrix(connect, 1, num, 1, 7);
    free_alc(num, nbond, AtomAlc, xyzAlc, ibase, 0, linkage);
}

static void alc2img_usage(void)
{
    help3dna_usage("alc2img");
}

static void set_defaults(struct_args * args)
{
    strcpy(args->inpfile, "");
    strcpy(args->outfile, "");
    strcpy(args->pdb, "");
    args->mol = FALSE;
    args->scale_factor = 0.0;
}

static void alc2img_cmdline(int argc, char *argv[], struct_args * args, long *opts)
{
    char *pdbdft = "ALC:A:1";
    long i, j;

    if (argc < 3)
        alc2img_usage();

    set_defaults(args);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        if (lux_ncmatch(argv[i], "^--?mo")) {
            if (strpbrk(argv[i], "ex3"))  /* v3000 -- *ex*tended molfile */
                args->mol = 3;
            else
                args->mol = 2;  /* default to v2000 */
            continue;
        }

        if (str_pmatch(argv[i], "-pdb")) {
            if (strchr(argv[i], '=') == NULL)
                strcpy(args->pdb, pdbdft);
            else {
                get_strvalue(argv[i], args->pdb, 0);
                if (is_empty_string(args->pdb))
                    strcpy(args->pdb, pdbdft);
            }
            continue;
        }

        upperstr(argv[i]);

        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'F')
                opts[0] = 1;
            else if (argv[i][j] == 'R')
                opts[0] = 2;
            else if (argv[i][j] == 'B')
                opts[1] = 1;
            else if (argv[i][j] == 'C')
                opts[2] = 1;
            else if (argv[i][j] == 'L')
                opts[3] = 1;
            else if (argv[i][j] == 'A')
                opts[5] = 1;
            else if (argv[i][j] == 'G')
                opts[6] = 1;
            else if (argv[i][j] == 'I')
                opts[7] = 1;
            else if (argv[i][j] == 'N')
                opts[8] = 1;
            else if (argv[i][j] == 'U')
                opts[9] = 1;
            else if (argv[i][j] == 'M')
                opts[10] = 1;
            else if (argv[i][j] == 'P')
                strcpy(args->pdb, pdbdft);
            else if (argv[i][j] == 'S') {
                if (sscanf(argv[i], "-S=%lf", &args->scale_factor) != 1) {
                    fprintf(stderr, "wrong format for setting scale factor\n");
                    alc2img_usage();
                } else
                    break;  /* -s=factor not combined with others */
            } else
                alc2img_usage();
    }
    if (argc == i + 2) {
        strcpy(args->inpfile, argv[i]);
        strcpy(args->outfile, argv[i + 1]);
    } else
        alc2img_usage();
    if (opts[10])  /* is 4 + minor/major faces, set -i option */
        opts[7] = 1;
    if (args->mol && !is_empty_string(args->pdb)) {
        fprintf(stderr, "Both -mol and -pdb options are set; ignore -pdb ...");
        strcpy(args->pdb, "");
    }
}

/* ALCHEMY schematic presentation to PS[0], XFIG[1], Raster3D[2] image */
int main(int argc, char *argv[])
{
    struct_args args;
    long opts[BUF512] = { 0 };  /* all initialized to zeros */
    /* 0: 0 for PS, 1 (F) for XFIG, and 2 (R) for Raster3D
     * 1: (B) frame_box
     * 2: (C) is_color
     * 3: (L) oline
     * 4: peptide (NOT used here, for compatibility with PDB2IMG)
     * 5: (A) helix_axis
     * 6: (G) ref_axis
     * 7: (I) same_faces
     * 8: (N) no_header
     * 9: (U) updown
     * 10: (M) isomM */

    set_my_globals(argv[0]);

    alc2img_cmdline(argc, argv, &args, opts);

    if (args.mol)
        cvt2mol(args.mol, args.inpfile, args.outfile);
    else if (!is_empty_string(args.pdb))
        cvt2pdb(&args);
    else
        process_alc(args.inpfile, args.outfile, args.scale_factor, opts);

    clear_my_globals();

    return 0;
}
