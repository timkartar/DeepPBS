#include "x3dna.h"

typedef struct {
    char parfile[BUF512];
    char inpfile[BUF512];
    char outfile[BUF512];
    long side_view;
    long helix_axis;
    long ref_axis;
    double axis_len;
    long is_helical;
    long n1;
    long n2;
    long bp_str;
    long frame;
    long reverse;  /* as in rebuild to align structure to a frame */
} struct_args;

static void write_frame_rotmat(long side_view, double *morg, double **mst)
{
    double **cmbrot;
    long i, j;

    cmbrot = dmatrix(1, 4, 1, 3);

    cpxyz(morg, cmbrot[4]);
    for (i = 1; i <= 3; i++)
        for (j = 1; j <= 3; j++)
            cmbrot[i][j] = mst[i][j];

    if (side_view)
        get_side_view(1, 3, cmbrot);

    write_rotmat(cmbrot);

    free_dmatrix(cmbrot, 1, 4, 1, 3);
}

static void frame_mol_usage(void)
{
    help3dna_usage("frame_mol");
}

static void set_defaults(struct_args * args)
{
    strcpy(args->parfile, "");
    strcpy(args->inpfile, "");
    strcpy(args->outfile, "");
    args->side_view = FALSE;
    args->helix_axis = FALSE;
    args->ref_axis = FALSE;
    args->axis_len = AXIS_LENGTH;
    args->is_helical = FALSE;
    args->n1 = 0;
    args->n2 = 0;
    args->bp_str = FALSE;
    args->frame = FALSE;
    args->reverse = FALSE;  /* rebuild, align to */
}

static void frame_mol_cmdline(int argc, char *argv[], struct_args * args)
{
    long i, j, n = 0;

    if (argc < 4)
        frame_mol_usage();

    set_defaults(args);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'M')
                args->side_view = TRUE;
            else if (argv[i][j] == 'A')
                args->helix_axis = TRUE;
            else if (argv[i][j] == 'G')
                args->ref_axis = TRUE;
            else if (argv[i][j] == 'L') {
                if (sscanf(argv[i], "-L=%lf", &args->axis_len) != 1) {
                    fprintf(stderr, "wrong format for setting axis length <%s>\n", argv[i]);
                    frame_mol_usage();
                } else
                    break;  /* -l=axis_length not combined with others */
            } else if (argv[i][j] == 'S')
                args->bp_str = TRUE;
            else if (argv[i][j] == 'F')
                args->frame = TRUE;
            else if (argv[i][j] == 'R')
                args->reverse = TRUE;
            else if (argv[i][j] == 'X')
                args->is_helical = TRUE;
            else if (isdigit((int) argv[i][j])) {
                n = sscanf(argv[i], "-%ld,%ld", &args->n1, &args->n2);
                if (n == 1) {
                    args->n2 = 0;  /* just to make sure */
                    if (args->is_helical) {
                        fprintf(stderr, "two numbers -n1,n2 required to"
                                " define a the helical frame\n");
                        frame_mol_usage();
                    }
                } else if (n == 2 || args->frame);  /* okay */
                else
                    frame_mol_usage();
                break;
            } else  /* wrong options or -h */
                frame_mol_usage();
    }
    if (argc == i + 3 && (n || args->frame)) {
        strcpy(args->parfile, argv[i]);
        strcpy(args->inpfile, argv[i + 1]);
        strcpy(args->outfile, argv[i + 2]);
    } else if (argc == i + 2 && args->bp_str) {
        strcpy(args->parfile, argv[i]);
        strcpy(args->outfile, argv[i + 1]);
    } else
        frame_mol_usage();
}

/* attach helix axes to ALCHEMY file */
static void attach_helix(long nbpm1, double **orien, double **org, long *num,
                         char **AtomName, long *ibase, double **xyz, long *nbond, long **linkage)
{
    double half_rise;
    double morg[4], pars[7];
    double **mst, **rot1, **rot2;
    long i, j, k;

    rot1 = dmatrix(1, 3, 1, 3);
    rot2 = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= nbpm1; i++) {
        orien2mst(orien[i], 0, rot1);
        orien2mst(orien[i + 1], 0, rot2);

        helical_par(rot1, org[i], rot2, org[i + 1], pars, mst, morg);
        k = *num + 2 * i - 1;
        half_rise = 0.5 * pars[3];
        for (j = 1; j <= 3; j++) {
            xyz[k][j] = morg[j] - half_rise * mst[j][3];
            xyz[k + 1][j] = morg[j] + half_rise * mst[j][3];
        }
        for (j = k; j <= k + 1; j++) {
            strcpy(AtomName[j], "OH");
            ibase[j] = NON_WC_IDX;  /* non-common base atoms */
        }
        linkage[*nbond + i][1] = k;
        linkage[*nbond + i][2] = k + 1;
    }
    *num += 2 * nbpm1;
    *nbond += nbpm1;

    free_dmatrix(rot1, 1, 3, 1, 3);
    free_dmatrix(rot2, 1, 3, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
}

/* generate base-pair block representation from reference frame */
static void ref_duplex_block(long num_bp, char **bp_seq, double **org, double **orien,
                             char *outfile)
{
    char bpi[3], BDIR[BUF512], **AtomName, **tAtomName;
    double **rmtx, **xyz, **txyz;
    long i, ik, ia = 0, ib = 0, j, k, bidx, nbond, num, tnbond, tnum;
    long *ibase, *tibase, **linkage, **tlinkage;

    get_BDIR(BDIR, "Block_BP.alc");
    strcat(BDIR, "Block_BP.alc");
    get_alc_nums(BDIR, &num, &nbond);

    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 1, 2);

    tnum = (num + 1) * num_bp;  /* including bp origin */
    tnbond = (nbond + 1) * num_bp - 1;  /* bp origin linkage */
    tAtomName = cmatrix(1, tnum, 0, 2);
    txyz = dmatrix(1, tnum, 1, 3);
    tibase = lvector(1, tnum);
    tlinkage = lmatrix(1, tnbond, 1, 2);

    rmtx = dmatrix(1, 3, 1, 3);

    read_alc(BDIR, &num, &nbond, AtomName, xyz, ibase, linkage);

    for (i = 1; i <= num_bp; i++) {
        orien2mst(orien[i], 0, rmtx);  /* get rotation matrix */

        sprintf(bpi, "%c%c", toupper((int) bp_seq[i][1]), toupper((int) bp_seq[i][2]));
        bidx = basepair_idx(bpi);

        for (j = 1; j <= num; j++) {
            ik = ia + j;
            strcpy(tAtomName[ik], AtomName[j]);
            tibase[ik] = bidx;
            for (k = 1; k <= 3; k++)
                txyz[ik][k] = dot(xyz[j], rmtx[k]) + org[i][k];
        }

        for (j = 1; j <= nbond; j++) {
            ik = ib + j;
            for (k = 1; k <= 2; k++)
                tlinkage[ik][k] = ia + linkage[j][k];
        }

        ia += num;
        ib += nbond;

    }

    cnct_org(num_bp, ia, ib, tAtomName, txyz, tibase, tlinkage, org);

    write_alc(tnum, tnbond, tAtomName, txyz, tibase, tlinkage, outfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_alc(tnum, tnbond, tAtomName, txyz, tibase, 1, tlinkage);
    free_dmatrix(rmtx, 1, 3, 1, 3);
}

/* undocumented: generate reference frame axes for illustration purpose */
static void pdbalc_frame(struct_args * args, long num_bp, double **org, double **orien)
{
    static char *ANPDB[4] = { " OO ", " OX ", " OY ", " OZ " };
    static char ANALC[3] = "OH";
    static char *ref_symbol[4] = { "OO", "OX", "OY", "OZ" };
    double d[4][3];
    static double ref_xyz[4][3] = {
        {0.0, 0.0, 0.0},  /* origin */
        {1.0, 0.0, 0.0},  /* x-axis */
        {0.0, 1.0, 0.0},  /* y-axis */
        {0.0, 0.0, 1.0}  /* z-axis */
    };
    long i, j, k, vnum0, vnum, nbond;
    long ipp = 0, ipa = 0, ipn = 0, ir = 0;
    FILE *fpdb, *falc;

    fpdb = open_file(args->inpfile, "w");
    falc = open_file(args->outfile, "w");

    vnum = vnum0 = 6 * num_bp;
    nbond = 3 * num_bp;
    if (args->ref_axis) {  /* with global reference frame */
        vnum += 4;
        nbond += 3;
    }

    /* Check later on with write_pdb/alc */
    fprintf(falc, "%5ld ATOMS, %5ld BONDS\n", vnum, nbond);

    for (i = 1; i <= num_bp; i++) {
        for (j = 1; j <= 3; j++) {
            d[0][j - 1] = org[i][j];
            for (k = 1; k <= 3; k++)
                d[j][k - 1] = BOND_UPPER_LIMIT * orien[i][(j - 1) * 3 + k] + org[i][k];
        }
        ir++;
        for (j = 0; j <= 3; j++)
            fprintf(fpdb, "ATOM  %5ld %4s %3s %c%4ld    %8.3f%8.3f%8.3f\n", ++ipp,
                    ANPDB[j], "AXE", 'Z', ir, d[j][0], d[j][1], d[j][2]);
        for (j = 1; j <= 3; j++) {
            fprintf(falc, "%5ld %2s   %9.4f%9.4f%9.4f\n", ++ipa, ANALC, d[0][0], d[0][1], d[0][2]);  /* origin */
            fprintf(falc, "%5ld %2s   %9.4f%9.4f%9.4f\n", ++ipa, ANALC, d[j][0], d[j][1], d[j][2]);  /* x, y, z tips */
        }
    }

    if (args->ref_axis) {  /* add global reference frame */
        ir++;
        for (j = 0; j <= 3; j++) {
            fprintf(fpdb, "ATOM  %5ld %4s %3s %c%4ld    ", ++ipp, ANPDB[j], "AXE", 'Z', ir);
            fprintf(falc, "%5ld %2s   ", ++ipa, ref_symbol[j]);
            for (k = 0; k <= 2; k++) {
                fprintf(fpdb, "%8.3f", args->axis_len * ref_xyz[j][k]);
                fprintf(falc, "%9.4f", args->axis_len * ref_xyz[j][k]);
            }
            fprintf(fpdb, "\n");
            fprintf(falc, "\n");
        }
    }

    for (i = 1; i <= ir; i++) {
        k = (i - 1) * 4;
        fprintf(fpdb, "CONECT%5ld%5ld%5ld%5ld\n", k + 1, k + 2, k + 3, k + 4);
        fprintf(fpdb, "CONECT%5ld%5ld\n", k + 2, k + 1);
        fprintf(fpdb, "CONECT%5ld%5ld\n", k + 3, k + 1);
        fprintf(fpdb, "CONECT%5ld%5ld\n", k + 4, k + 1);
    }

    ipn = 0;
    for (i = 1; i <= vnum0; i += 2)
        fprintf(falc, "%5ld %5ld %5ld\n", ++ipn, i, i + 1);
    if (args->ref_axis)
        for (j = 1; j <= 3; j++)
            fprintf(falc, "%5ld %5ld %5ld\n", ++ipn, vnum0 + 1, vnum0 + j + 1);

    close_file(fpdb);
    close_file(falc);
}

/* get mst and org */
static void get_mstorg(struct_args * args, double **org, double **orien, double *morg,
                       double **mst)
{
    double pars[7];
    double **rot1, **rot2;
    long i, bz = 1, bz_junction = 0, z_step = 0;

    /* define the new reference frame  */
    rot1 = dmatrix(1, 3, 1, 3);
    rot2 = dmatrix(1, 3, 1, 3);
    if (!args->n2) {  /* single base-pair */
        cpxyz(org[args->n1], morg);
        orien2mst(orien[args->n1], 0, mst);
    } else {  /* two base-pairs */
        orien2mst(orien[args->n1], 0, rot1);
        orien2mst(orien[args->n2], 0, rot2);

        bz_check(rot1, org[args->n1], rot2, org[args->n2], bz, &bz_junction, &z_step);
        if (args->is_helical)
            helical_par(rot1, org[args->n1], rot2, org[args->n2], pars, mst, morg);
        else
            bpstep_par(rot1, org[args->n1], rot2, org[args->n2], pars, mst, morg);
        fprintf(stderr, "the six parameters are:\n");
        for (i = 1; i <= 6; i++)
            fprintf(stderr, "%8.2f", pars[i]);
        fprintf(stderr, "\n");
    }
    free_dmatrix(rot1, 1, 3, 1, 3);
    free_dmatrix(rot2, 1, 3, 1, 3);
}

static void frame_pdbfile(long num, double *morg, double **mst, struct_args * args)
{
    char *ChainID, **AtomName, **ResName, **Miscs;
    double **xyz;
    long *ResSeq;

    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);

    Gvars.AtomName0 = cmatrix(1, num, 0, 4);
    Gvars.ResName0 = cmatrix(1, num, 0, 3);
    Gvars.Name0 = TRUE;

    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);

    read_pdb(args->inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");
    Gvars.Name0 = FALSE;

    if (args->reverse)
        frame_xyz(args->side_view, morg, mst, num, xyz);
    else
        change_xyz(args->side_view, morg, mst, num, xyz);
    write_pdb(num, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, args->outfile);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);

    free_cmatrix(Gvars.AtomName0, 1, DUMMY, 0, DUMMY);
    free_cmatrix(Gvars.ResName0, 1, DUMMY, 0, DUMMY);
}

static void frame_alcfile(long num_bp, double **org, double **orien, double *morg,
                          double **mst, struct_args * args)
{
    char **AtomName;
    double **xyz;
    long nbond, nbpm1, num, tnbond, tnum;
    long *ibase, **linkage;

    nbpm1 = num_bp - 1;

    get_alc_nums(args->inpfile, &num, &nbond);
    tnum = num + 2 * nbpm1 + 4;
    tnbond = nbond + nbpm1 + 3;

    AtomName = cmatrix(1, tnum, 0, 2);
    xyz = dmatrix(1, tnum, 1, 3);
    ibase = lvector(1, tnum);
    linkage = lmatrix(1, tnbond, 1, 2);

    read_alc(args->inpfile, &num, &nbond, AtomName, xyz, ibase, linkage);

    if (args->helix_axis)
        attach_helix(nbpm1, orien, org, &num, AtomName, ibase, xyz, &nbond, linkage);

    if (args->reverse)
        frame_xyz(args->side_view, morg, mst, num, xyz);
    else
        change_xyz(args->side_view, morg, mst, num, xyz);

    if (args->ref_axis)
        add_3axes(&num, AtomName, ibase, xyz, &nbond, linkage, args->side_view, args->axis_len);

    write_alc(num, nbond, AtomName, xyz, ibase, linkage, args->outfile);

    free_alc(tnum, tnbond, AtomName, xyz, ibase, 1, linkage);
}

typedef struct {
    long ds;
    long num_bp;
    char **bp_seq;
    double **org;
    double **orien;
} struct_frames;

static void init_frames(struct_frames * frames)
{
    frames->ds = FALSE;
    frames->num_bp = FALSE;
    frames->bp_seq = NULL;
    frames->org = NULL;
    frames->orien = NULL;
}

static void free_frames(struct_frames * frames)
{
    if (frames->bp_seq)
        free_cmatrix(frames->bp_seq, 1, DUMMY, 1, DUMMY);

    if (frames->org)
        free_dmatrix(frames->org, 1, DUMMY, 1, DUMMY);

    if (frames->orien)
        free_dmatrix(frames->orien, 1, DUMMY, 1, DUMMY);
}

/* read in reference frames: origin and orientation */
static void read_frames(char *parfile, struct_frames * frames)
{
    char *format = "%lf %lf %lf", stmp[BUF512];
    char str[BUF512];
    long i, num1 = 0, num2 = 0;
    FILE *fp;

    fp = open_file(parfile, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", &frames->num_bp) != 1)
        fatal("error reading ref-frame file: number of base-pairs\n");
    if (frames->num_bp <= 0)
        fatal("number of base-pairs [%ld] <= 0\n", frames->num_bp);

    frames->org = dmatrix(1, frames->num_bp, 1, 3);
    frames->orien = dmatrix(1, frames->num_bp, 1, 9);
    frames->bp_seq = cmatrix(1, frames->num_bp, 1, 2);

    for (i = 1; i <= frames->num_bp; i++) {
        if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%*s %*d %s", stmp) != 1)
            fatal("error reading title line\n");
        else {
            frames->bp_seq[i][1] = stmp[0];
            if (strlen(stmp) == 3) {
                num2++;
                frames->bp_seq[i][2] = stmp[2];
            } else
                num1++;
        }
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &frames->org[i][1], &frames->org[i][2], &frames->org[i][3]) != 3)
            fatal("error reading origin line\n");
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &frames->orien[i][1], &frames->orien[i][2],
                   &frames->orien[i][3]) != 3)
            fatal("error reading x direction cosines\n");
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &frames->orien[i][4], &frames->orien[i][5],
                   &frames->orien[i][6]) != 3)
            fatal("error reading y direction cosines\n");
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &frames->orien[i][7], &frames->orien[i][8],
                   &frames->orien[i][9]) != 3)
            fatal("error reading z direction cosines\n");
    }
    if (num1 && num2)
        fatal("duplex and single helix mixed\n");
    if ((num1 && num1 != frames->num_bp) || (num2 && num2 != frames->num_bp))
        fatal("reference frame data file not consistent\n");
    frames->ds = num2 ? 2 : 1;

    close_file(fp);
}

/* A 3DNA generated structure is always oriented w.r.t. the first
   base-pair (or base). A data file containing the reference frame
   (default name "ref_frames.dat") for each base-pair (or base) is
   also available, which is used here to set the new reference frame
   for the structure. The new frame could that of any base-pair (or
   base) or the middle frame or middle helical frame define by any two
   base-pairs (or bases) */
int main(int argc, char *argv[])
{
    double morg[4], **mst;
    long num, num_bp;
    struct_args args;
    struct_frames frames;

    set_my_globals(argv[0]);

    frame_mol_cmdline(argc, argv, &args);

    init_frames(&frames);
    read_frames(args.parfile, &frames);
    num_bp = frames.num_bp;

    if (args.bp_str && frames.ds == 2) {
        ref_duplex_block(num_bp, frames.bp_seq, frames.org, frames.orien, args.outfile);
        return 0;
    } else if (args.frame) {
        pdbalc_frame(&args, num_bp, frames.org, frames.orien);
        return 0;
    }

    if (!lval_in_range(args.n1, 1, num_bp))
        fatal("base-pair number out of range\n");
    if (args.n2 && (!lval_in_range(args.n2, 1, num_bp) || args.n2 <= args.n1))
        fatal("base-pair number out of range\n");

    /* define the new reference frame  */
    mst = dmatrix(1, 3, 1, 3);
    get_mstorg(&args, frames.org, frames.orien, morg, mst);
    write_frame_rotmat(args.side_view, morg, mst);

    num = number_of_atoms(args.inpfile, 1, "*");
    if (num)  /* is a PDB file */
        frame_pdbfile(num, morg, mst, &args);
    else
        frame_alcfile(num_bp, frames.org, frames.orien, morg, mst, &args);

    free_frames(&frames);
    free_dmatrix(mst, 1, 3, 1, 3);

    clear_my_globals();

    return 0;
}
