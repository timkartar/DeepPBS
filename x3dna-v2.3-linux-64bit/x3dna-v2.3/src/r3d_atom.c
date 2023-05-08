#include "x3dna.h"

typedef struct {
    char hbfile[BUF512];
    char lkgfile[BUF512];
    char inpfile[BUF512];
    char outfile[BUF512];
    double rodrad;
    double ballrad;
    double scale_factor;
} struct_args;

static void r3d_atom_usage(void)
{
    help3dna_usage("r3d_atom");
}

static void set_defaults(struct_args * args)
{
    strcpy(args->hbfile, "");
    strcpy(args->lkgfile, "");
    strcpy(args->inpfile, "");
    strcpy(args->outfile, "");
    args->rodrad = 0.2;
    args->ballrad = 0.0;
    args->scale_factor = 0.0;
}

static void r3d_atom_cmdline(int argc, char *argv[], struct_args * args, long *opts)
{
    char str[BUF512];  /* keep a copy of the original */
    long i, j;

    if (argc < 3)
        r3d_atom_usage();

    set_defaults(args);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        strcpy(str, argv[i]);
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'N')
                opts[1] = 1;
            else if (argv[i][j] == 'O')
                opts[2] = 1;
            else if (argv[i][j] == 'D')
                opts[3] = 1;
            else if (argv[i][j] == 'C')
                opts[4] = 1;
            else if (argv[i][j] == 'P')
                opts[5] = 1;
            else if (argv[i][j] == 'G')
                opts[6] = 1;
            else if (argv[i][j] == 'L')
                opts[7] = 1;
            else if (argv[i][j] == 'Z')
                opts[8] = 1;
            else if (argv[i][j] == 'X') {
                if (sscanf(str, "-X=%s", args->hbfile) != 1 &&
                    sscanf(str, "-x=%s", args->hbfile) != 1) {
                    fprintf(stderr, "wrong format for reading H-bond list\n");
                    r3d_atom_usage();
                } else {
                    opts[3] = 1;
                    break;  /* -x=hbond not combined with others */
                }
            } else if (argv[i][j] == 'Y') {
                if (sscanf(str, "-Y=%s", args->lkgfile) != 1 &&
                    sscanf(str, "-y=%s", args->lkgfile) != 1) {
                    fprintf(stderr, "wrong format for reading bond list\n");
                    r3d_atom_usage();
                } else
                    break;  /* -y=bond list file not combined with others */
            } else if (argv[i][j] == 'S') {
                if (sscanf(argv[i], "-S=%lf", &args->scale_factor) != 1) {
                    fprintf(stderr, "wrong format for setting scale factor\n");
                    r3d_atom_usage();
                } else
                    break;
            } else if (argv[i][j] == 'R') {
                if (sscanf(argv[i], "-R=%lf", &args->rodrad) != 1) {
                    fprintf(stderr, "wrong format for cylinder radius\n");
                    r3d_atom_usage();
                } else
                    break;
            } else if (argv[i][j] == 'B') {
                if (sscanf(argv[i], "-B=%lf", &args->ballrad) != 1) {
                    fprintf(stderr, "wrong format for sphere radius\n");
                    r3d_atom_usage();
                } else
                    break;
            } else
                r3d_atom_usage();
    }
    if (argc == i + 2) {
        strcpy(args->inpfile, argv[i]);
        strcpy(args->outfile, argv[i + 1]);
    } else
        r3d_atom_usage();
    if (args->rodrad <= 0.0 && args->ballrad <= 0.0 && !opts[2])
        args->rodrad = 0.20;
}

/* stick presentation materials to be rendered by Raster3d */
static void stick_model(long nobase, char *lkgfile, long num, long *RY, long *idx,
                        char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                        char **Miscs, double **xyz, double rodrad, long *colidx,
                        double **colrgb, FILE * fp)
{
    long i, ia, ib, nbond = 0, nbond_estimated, bond_lmt;
    long *baseatom, **linkage;
    double dxyz[4];

    if (rodrad <= 0.0)
        return;

    bond_lmt = lround(NBOND_FNUM * num);
    nbond_estimated = 5 * bond_lmt;  /* to account for overlaps */
    linkage = lmatrix(1, nbond_estimated, 1, 2);

    fprintf(fp, "###\n### The following section is for stick/rod model\n");
    if (!strcmp(lkgfile, "")) {  /* no bond linkage file supplied */
        atom_linkage(1, num, idx, xyz, Miscs, ChainID, nbond_estimated, &nbond, linkage);
        write_pdb(num, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, SNUM_FILE);
        write_lkglist(nbond, linkage, AtomName, ResName, ChainID, ResSeq, Miscs);
    } else
        read_lkginfo(lkgfile, num, &nbond, linkage);

    if (nbond > bond_lmt)
        fprintf(stderr, "too many linkages: *** possibly overlap ***\n");

    baseatom = lvector(1, num);
    if (rodrad > 0) {  /* in case no bonds are needed */
        if (nobase)  /* do not draw bonds within base atoms */
            for (i = 1; i <= num; i++) {
                if (RY[i] < 0)
                    continue;
                baseatom[i] = is_baseatom(AtomName[i]);
            }
        for (i = 1; i <= nbond; i++) {
            ia = linkage[i][1];
            ib = linkage[i][2];
            if (nobase && baseatom[ia] && baseatom[ib])
                continue;
            if (colidx[ia] == colidx[ib])  /* same kind of colors */
                r3d_rod(3, xyz[ia], xyz[ib], rodrad, colrgb[ia], fp);
            else {
                avexyz(xyz[ia], xyz[ib], dxyz);  /* middle point of atoms ia & ib */
                r3d_rod(3, xyz[ia], dxyz, rodrad, colrgb[ia], fp);
                r3d_rod(3, dxyz, xyz[ib], rodrad, colrgb[ib], fp);
            }
        }
    }

    write_alc(num, nbond, AtomName, xyz, NULL, linkage, ATOMALC_FILE);

    free_lmatrix(linkage, 1, nbond_estimated, 1, 2);
    free_lvector(baseatom, 1, num);
}

/* reset color for gray image */
static void gray_image(long is_gray, double **atom_col, double **base_col, double *hb_col)
{
    long i;

    if (is_gray) {  /* gray image */
        cpxyz(base_col[NBASECOL], hb_col);
        for (i = 0; i <= NATOMCOL; i++)
            cpxyz(base_col[NBASECOL], atom_col[i]);
        for (i = 0; i <= NBASECOL; i++)
            cpxyz(base_col[NBASECOL], base_col[i]);
    }
}

/* color by base residue */
static void col_residue(long num, long num_residue, long *RY, long **seidx, long *idx,
                        double **atom_col, double **base_col, long *ibase,
                        long by_residue, long *colidx, double **colrgb)
{
    long i, j;

    for (i = 1; i <= num; i++) {
        colidx[i] = idx[i];
        cpxyz(atom_col[idx[i]], colrgb[i]);
    }

    if (by_residue)  /* colored by base residues or by atoms */
        for (i = 1; i <= num_residue; i++)
            if (RY[i] >= 0) {  /* base residue */
                for (j = seidx[i][1]; j <= seidx[i][2]; j++) {
                    colidx[j] = ibase[i];
                    cpxyz(base_col[ibase[i]], colrgb[j]);
                }
            }
}

/* render H-bonding information */
static void render_hbonds(long num, long num_residue, char *bseq, long **seidx, long *idx,
                          long *RY, char **AtomName, char **ResName, char *ChainID,
                          long *ResSeq, char **Miscs, double **xyz, long pwise,
                          double *width3, double rodrad, char *hbfile, double *hb_col, FILE * fp)
{
    double hb_width;
    long i, num_hbond = 0, **hb_linkage;

    hb_width = (rodrad <= 0.0) ? width3[1] : rodrad;
    hb_linkage = lmatrix(1, num, 1, 3);
    if (!strcmp(hbfile, ""))  /* no hbfile inputed */
        hbond_pdb(num, num_residue, bseq, seidx, idx, RY, AtomName, ResName, ChainID,
                  ResSeq, Miscs, xyz, &num_hbond, hb_linkage, pwise);
    else
        read_hbinfo(hbfile, num, &num_hbond, hb_linkage);

    fprintf(fp, "###\n### The following section is for H-bonds\n");
    for (i = 1; i <= num_hbond; i++)
        r3d_dash(xyz[hb_linkage[i][2]], xyz[hb_linkage[i][3]], hb_width, hb_col, fp);

    free_lmatrix(hb_linkage, 1, num, 1, 3);
}

static void r3d_input(struct_args * args, long *opts)
{
    char label_style[BUF512], *ChainID, *bseq, **AtomName, **ResName, **Miscs;
    double width3[4], hb_col[5];
    double **atom_col, **base_col, **colrgb, **xyz;
    long num, num_ring, num_residue;
    long *ibase, *idx, *colidx, *ResSeq, *RY, **ring_atom, **seidx;
    FILE *fp;

    num = number_of_atoms(args->inpfile, 1, "*");
    if (!num)
        fatal("PDB file without ATOM/HETATM records\n");
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(args->inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");

    /* get the numbering information of each residue */
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

    bseq = cvector(1, num_residue);
    RY = lvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);

    ring_atom = lmatrix(1, num_residue, 1, 10);
    all_bring_atoms(num_residue, RY, seidx, AtomName, &num_ring, ring_atom);

    ibase = lvector(1, num_residue);
    base_idx(num_residue, bseq, ibase, 0);
    idx = lvector(1, num);
    atom_idx(num, AtomName, Miscs, idx);

    atom_col = dmatrix(0, NATOMCOL, 1, 3);
    base_col = dmatrix(0, NBASECOL, 1, 3);
    get_r3dpars(base_col, hb_col, width3, atom_col, label_style);

    gray_image(opts[6], atom_col, base_col, hb_col);

    colidx = lvector(1, num);
    colrgb = dmatrix(1, num, 1, 3);
    col_residue(num, num_residue, RY, seidx, idx, atom_col, base_col, ibase,
                opts[4], colidx, colrgb);

    fp = open_file(args->outfile, "w");
    raster3d_header(num, xyz, args->scale_factor, opts[1], 0, fp);

    cpk_model(num, idx, xyz, args->ballrad, colrgb, fp);
    stick_model(opts[8], args->lkgfile, num, RY, idx, AtomName, ResName, ChainID, ResSeq,
                Miscs, xyz, args->rodrad, colidx, colrgb, fp);

    if (opts[2])
        fill_base_ring(num_residue, num_ring, ring_atom, xyz, ibase, bseq,
                       base_col, label_style, opts[7], ResSeq, fp);

    if (opts[3] || opts[5])  /* draw H-bonds between bases */
        render_hbonds(num, num_residue, bseq, seidx, idx, RY, AtomName, ResName, ChainID,
                      ResSeq, Miscs, xyz, opts[5], width3, args->rodrad, args->hbfile, hb_col, fp);

    close_file(fp);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cvector(bseq, 1, num_residue);
    free_lvector(RY, 1, num_residue);
    free_lmatrix(ring_atom, 1, num_residue, 1, 10);
    free_lvector(ibase, 1, num_residue);
    free_lvector(idx, 1, num);
    free_dmatrix(atom_col, 0, NATOMCOL, 1, 3);
    free_dmatrix(base_col, 0, NBASECOL, 1, 3);
    free_lvector(colidx, 1, num);
    free_dmatrix(colrgb, 1, num, 1, 3);
}

int main(int argc, char *argv[])
{
    struct_args args;
    long opts[BUF512] = { 0 };
    /* 1: (N) no-header
     * 2: (O) fill base ring
     * 3: (D) draw-hbonds;
     * 4: (C) color-residue
     * 5: (P) pairwise when finding hbonds (multiplets from find_pair)
     * 6: (G) gray image
     * 7: (L) label base residue in center of 6-member ring
     * 8: (Z) do not draw bonds within base atoms for "blocview"
     * -x=hbond_list
     * -y=bond_list */

    set_my_globals(argv[0]);

    r3d_atom_cmdline(argc, argv, &args, opts);

    if (strcmp(args.hbfile, HB_FILE))
        remove_file(HB_FILE);
    if (strcmp(args.lkgfile, LKG_FILE))
        remove_file(LKG_FILE);
    remove_file(ATOMALC_FILE);
    remove_file(SNUM_FILE);

    r3d_input(&args, opts);

    clear_my_globals();

    return 0;
}
