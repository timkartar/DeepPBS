#include "x3dna.h"

typedef struct {
    char rotfile[BUF512];
    char inpfile[BUF512];
    char outfile[BUF512];
    long fnum;
    long to_center;
    long pmi;
} struct_args;

static void rotmol_usage(void)
{
    help3dna_usage("rotate_mol");
}

static void set_defaults(struct_args * args)
{
    strcpy(args->rotfile, "");
    strcpy(args->inpfile, "");
    strcpy(args->outfile, "");
    args->fnum = 0;
    args->to_center = FALSE;
    args->pmi = 2;  /* pmi based on nucleotides */
}

static void rotate_mol_cmdline(long argc, char *argv[], struct_args * args)
{
    char str[BUF512];  /* keep a copy of the original */
    long i, j;

    if (argc < 3)
        rotmol_usage();

    set_defaults(args);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        strcpy(str, argv[i]);
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++) {  /* skip - */
            if (argv[i][j] == 'C')
                args->to_center = TRUE;
            else if (argv[i][j] == 'A')
                args->pmi = 1;  /* all atoms */
            else if (argv[i][j] == 'N')
                args->pmi = 2;  /* nucleotides */
            else if (argv[i][j] == 'P')
                args->pmi = 3;  /* amino-acids */
            else if (argv[i][j] == 'B')
                args->pmi = 4;  /* bases */
            else if (argv[i][j] == 'R') {
                if (sscanf(str, "-R=%s", args->rotfile) != 1 &&
                    sscanf(str, "-r=%s", args->rotfile) != 1) {
                    fprintf(stderr, "wrong format for reading rotation angles\n");
                    rotmol_usage();
                } else {
                    args->fnum = 1;
                    break;  /* -r=rotfile not combined with others */
                }
            } else if (argv[i][j] == 'T') {
                if (sscanf(str, "-T=%s", args->rotfile) != 1 &&
                    sscanf(str, "-t=%s", args->rotfile) != 1) {
                    fprintf(stderr, "wrong format for reading transformation matrix\n");
                    rotmol_usage();
                } else {
                    args->fnum = 2;
                    break;  /* -t=transMtx file not combined with others */
                }
            } else
                rotmol_usage();
        }
    }
    if (argc == i + 2) {
        strcpy(args->inpfile, argv[i++]);
        strcpy(args->outfile, argv[i]);
    } else
        rotmol_usage();
}

/* get rotation matrix based on MolScript input (as from RasMol) */
static void molscr_rotmat(char *rotfile, double **rotmat)
{
    char axis, str[BUF512], str_by[BUF512], str_rotation[BUF512];
    double rot_ang, **roti, **temp;
    long naxis = 0;
    FILE *fp;

    fp = open_file(rotfile, "r");
    roti = dmatrix(1, 3, 1, 3);
    temp = dmatrix(1, 3, 1, 3);
    while (fgets(str, sizeof str, fp) != NULL) {
        upperstr(str);  /* change to UPPER case */
        if (sscanf(str, "%s %s", str_by, str_rotation) == 2 &&
            !strcmp(str_by, "BY") && !strcmp(str_rotation, "ROTATION")) {
            if (sscanf(str, "%*s %*s %c %lf", &axis, &rot_ang) != 2)
                fatal("wrong <by rotation> format\n");
            if (axis == 'X' || axis == 'Y' || axis == 'Z')
                naxis++;
            else
                fatal("none of x-, y- or z-axis\n");
            if (axis == 'X')
                rotx(-rot_ang, roti);
            else if (axis == 'Y')
                roty(-rot_ang, roti);
            else  /* axis == 'Z' */
                rotz(-rot_ang, roti);
            multi_matrix(rotmat, 3, 3, roti, 3, 3, temp);
            copy_dmatrix(temp, 3, 3, rotmat);
        }
    }
    close_file(fp);
    if (!naxis)
        fprintf(stderr, "NO rotation applied to the structure\n");

    free_dmatrix(roti, 1, 3, 1, 3);
    free_dmatrix(temp, 1, 3, 1, 3);
}

static void get_selected_xyz(long *num_sel, long j, double **xyz, double **sel_xyz,
                             char **AtomName, char **ResName, char *ChainID, long *ResSeq)
{
    (*num_sel)++;
    cpxyz(xyz[j], sel_xyz[*num_sel]);

    if (Gvars.VERBOSE)
        fprintf(stderr, "%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f\n", *num_sel,
                AtomName[j], ' ', ResName[j], ChainID[j], ResSeq[j], ' ',
                xyz[j][1], xyz[j][2], xyz[j][3]);
}

static void print_verbose_x_y_z(double *x, double *y, double *z)
{
    if (Gvars.VERBOSE) {
        fprintf(stderr, "x-axis: %10.6f %10.6f %10.6f\n", x[1], x[2], x[3]);
        fprintf(stderr, "y-axis: %10.6f %10.6f %10.6f\n", y[1], y[2], y[3]);
        fprintf(stderr, "z-axis: %10.6f %10.6f %10.6f\n", z[1], z[2], z[3]);
    }
}

static void flip_x_or_y(double *x, double *y, double **rotmat)
{
    double **tmpmtx, **rot180;

    tmpmtx = dmatrix(1, 3, 1, 3);
    copy_dmatrix(rotmat, 3, 3, tmpmtx);

    rot180 = dmatrix(1, 3, 1, 3);

    if (((fabs(x[2]) > fabs(x[1])) && x[2] < 0.0) ||  /* case #1 */
        ((fabs(y[2]) > fabs(y[1])) && y[2] > 0.0)) {  /* case #2 */
        rotx(180.0, rot180);
        multi_matrix(tmpmtx, 3, 3, rot180, 3, 3, rotmat);
        copy_dmatrix(rotmat, 3, 3, tmpmtx);  /* copy again for later reference */
        if (Gvars.VERBOSE)
            fprintf(stderr, "flipped vertically: rotx_180\n");
    }

    if (((fabs(y[1]) > fabs(y[2])) && y[1] > 0.0) ||  /* case #1 */
        ((fabs(x[1]) > fabs(x[2])) && x[1] < 0.0)) {  /* case #2 */
        roty(180.0, rot180);
        multi_matrix(tmpmtx, 3, 3, rot180, 3, 3, rotmat);
        if (Gvars.VERBOSE)
            fprintf(stderr, "flipped horizontally: roty_180\n");
    }

    free_dmatrix(tmpmtx, 1, 3, 1, 3);
    free_dmatrix(rot180, 1, 3, 1, 3);
}

/* double check base 1 frame to avoid ambiguity in overall orientation:
 * x-axis: pointing upwards; y-axis: pointing left; z-axis: point out */
static void check_base1(long num_residue, long *RY, double **orien, double **rotmat)
{
    double *x, *y, *z;
    long i;

    for (i = 1; i <= num_residue; i++)
        if (RY[i] >= 0)
            break;
    if (i > num_residue)  /* no nucleotides */
        return;

    x = &orien[i][0];
    y = &orien[i][3];
    z = &orien[i][6];
    print_verbose_x_y_z(x, y, z);

/*
see user_cases / pascal_rotate_mol for examples:
   ga_1.pdb -- needs rotx180, then roty180
           x-axis:     0.00    -1.00    -0.03
           y-axis:     1.00     0.00    -0.01
           z-axis:     0.01    -0.03     1.00

   ga_2.pdb -- desired orientation
           x-axis:     0.07     1.00     0.00
           y-axis:    -1.00     0.07    -0.03
           z-axis:    -0.03    -0.00     1.00
*/

    flip_x_or_y(x, y, rotmat);
}

/* double check the frame defined by the first three atoms to avoid
 * ambiguity in overall orientation:  x-axis: pointing upwards;
 * y-axis:  pointing left; z-axis: point out */
static void check_first_3atoms(long num, double **xyz, double **rotmat)
{
    double x[4], y[4], z[4], **refmat;

    refmat = dmatrix(1, 4, 1, 3);
    identity_matrix(refmat, 3);
    define_frame_by_3atoms(num, xyz, refmat);

    mtx_2_x_y_z(refmat, x, y, z);
    print_verbose_x_y_z(x, y, z);

    flip_x_or_y(x, y, rotmat);

    free_dmatrix(refmat, 1, 4, 1, 3);
}

/* get the new view of PMI */
static void pmi_view(double **rotmat, double **tmprot, long c1, long c2, long c3)
{
    long i;

    for (i = 1; i <= 3; i++) {
        tmprot[i][1] = rotmat[i][c1];  /* x-axis */
        tmprot[i][2] = rotmat[i][c2];  /* y-axis */
        tmprot[i][3] = rotmat[i][c3];  /* z-axis */
    }
}

/* get the principle moment of inertia. By default, x-axis points
   along the normal to the best plane (least variance), z-axis along
   the best line (most variance), and y-axis in between. Reorder it as
   [2 3 1] to get the most-extended view and make sure a(?)-axis has
   positive components: 1ber */
static void pmi_rotmat(long num, double **xyz, double **rotmat)
{
    double tx[4], ty[4], tz[4], **cov_mtx, **rot180;
    long i;

    if (!num)
        return;
    cov_mtx = dmatrix(1, 3, 1, 3);
    cov_matrix(xyz, xyz, num, 3, cov_mtx);
    jacobi(cov_mtx, 3, tx, rotmat);

    /* Make sure (:, 3) == cross((:, 1),  (:, 2)) */
    for (i = 1; i <= 3; i++) {
        tx[i] = rotmat[i][1];
        ty[i] = rotmat[i][2];
    }
    cross(tx, ty, tz);
    for (i = 1; i <= 3; i++)
        rotmat[i][3] = tz[i];

    pmi_view(rotmat, cov_mtx, 2, 3, 1);  /* 2-3-1 view */
    copy_dmatrix(cov_mtx, 3, 3, rotmat);

    rot180 = dmatrix(1, 3, 1, 3);
    if (cov_mtx[1][1] < 0.0) {  /* x-axis reversed */
        roty(180.0, rot180);
        multi_matrix(cov_mtx, 3, 3, rot180, 3, 3, rotmat);
        copy_dmatrix(rotmat, 3, 3, cov_mtx);  /* copy back */
    }
    if (cov_mtx[2][2] < 0.0) {  /* y-axis reversed */
        rotx(180.0, rot180);
        multi_matrix(cov_mtx, 3, 3, rot180, 3, 3, rotmat);
    }

    free_dmatrix(cov_mtx, 1, 3, 1, 3);
    free_dmatrix(rot180, 1, 3, 1, 3);
}

/* reorient the structure with regard to the first base reference frame if any.
   otherwise, use the first three atoms to define the reference frame */
static void reorient_xyz(long num, double **xyz, double *orien1, double *org1, double **refmat)
{
    if (orien1 == NULL)
        define_frame_by_3atoms(num, xyz, refmat);
    else {  /* with base reference */
        orien2mst(orien1, 0, refmat);
        cpxyz(org1, refmat[4]);
    }

    /* ignore translations */
    init_dvector(refmat[4], 1, 3, 0.0);

    change_xyz(0, refmat[4], refmat, num, xyz);
}

/* get rotation matrix based on PMI of a PDB file */
static void pdb_pmi(long num, char **AtomName, char **ResName, char *ChainID,
                    long *ResSeq, double **xyz, char **Miscs, long pmi, double *ave_xyz,
                    double **refmat, double **rotmat)
{
    char BDIR[BUF512], *bseq;
    double **orien, **org, **sel_xyz;
    long i, j, nt = 0, num_residue, num_sel = 0, *RY, **seidx;

    if (num < 3)
        return;

    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);
    bseq = cvector(1, num_residue);
    RY = lvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);

    for (i = 1; i <= num_residue; i++)
        if (RY[i] >= 0)
            nt++;

    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 3);
    get_BDIR(BDIR, "Atomic_A.pdb");
    base_frame(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID, ResSeq, Miscs,
               xyz, BDIR, orien, org);

    /* xyz coordinates are transformed via refmat; no translation */
    for (i = 1; i <= num_residue; i++)
        if (RY[i] >= 0) {
            reorient_xyz(num, xyz, orien[i], org[i], refmat);
            break;
        }
    if (i > num_residue)  /* no nucleotide */
        reorient_xyz(num, xyz, NULL, NULL, refmat);

    sel_xyz = dmatrix(1, num, 1, 3);
    for (i = 1; i <= num_residue; i++) {
        if ((pmi == 2 && RY[i] >= 0) ||  /* base residue */
            (pmi == 3 && RY[i] == -1)) {  /* protein */
            for (j = seidx[i][1]; j <= seidx[i][2]; j++)
                get_selected_xyz(&num_sel, j, xyz, sel_xyz, AtomName, ResName, ChainID, ResSeq);
        } else if (pmi == 4 && RY[i] >= 0) {
            for (j = seidx[i][1]; j <= seidx[i][2]; j++) {
                if (is_baseatom(AtomName[j]))
                    get_selected_xyz(&num_sel, j, xyz, sel_xyz, AtomName, ResName,
                                     ChainID, ResSeq);
            }
        }
    }

    if (pmi >= 2 && num_sel) {  /* with selected atoms */
        pmi_rotmat(num_sel, sel_xyz, rotmat);
        ave_dmatrix(sel_xyz, num_sel, 3, ave_xyz);
    } else {  /* all atoms used */
        pmi_rotmat(num, xyz, rotmat);
        ave_dmatrix(xyz, num, 3, ave_xyz);
    }

    /* sel_xyz[][] renewed to re-calculate frames */
    multi_matrix(xyz, num, 3, rotmat, 3, 3, sel_xyz);
    if (nt > 0) {  /* with nucleotide */
        base_frame(num_residue, bseq, seidx, RY, AtomName, ResName,
                   ChainID, ResSeq, Miscs, sel_xyz, BDIR, orien, org);
        check_base1(num_residue, RY, orien, rotmat);
    } else
        check_first_3atoms(num, sel_xyz, rotmat);

    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cvector(bseq, 1, num_residue);
    free_lvector(RY, 1, num_residue);
    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 3);
    free_dmatrix(sel_xyz, 1, num, 1, 3);
}

/* combine two transformation matrices: applied to PMI view
 * [(xyz - o1) * r1 - o2] * r2
 *   ====> [(xyz - o1) * r1 - o2] * (inv(r1) * r1 * r2)
 *   ====> [xyz - (o1 + o2 * inv(r1))] * (r1 * r2)
 *  O = o1 + o2 * inv(r1); R = r1 * r2
 * inv(r1) = transpose of r1
 * here r1: refmat; r2: rotmat */
static void cmb_transformation(double **refmat, double **rotmat)
{
    double **cmbrot;

    cmbrot = dmatrix(1, 4, 1, 3);

    multi_vec_Tmatrix(rotmat[4], 3, refmat, 3, 3, cmbrot[4]);
    sumxyz(refmat[4], cmbrot[4], cmbrot[4]);
    multi_matrix(refmat, 3, 3, rotmat, 3, 3, cmbrot);
    write_rotmat(cmbrot);

    free_dmatrix(cmbrot, 1, 4, 1, 3);
}

static void rot_pdbfile(long num, double **rotmat, struct_args * args)
{
    char *ChainID, **AtomName, **ResName, **Miscs;
    double **nxyz, **xyz, **refmat, **tmprot, ave_xyz[4];
    long *ResSeq;

    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);

    Gvars.AtomName0 = cmatrix(1, num, 0, 4);
    Gvars.ResName0 = cmatrix(1, num, 0, 3);
    Gvars.Name0 = TRUE;

    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    nxyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);

    read_pdb(args->inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");
    Gvars.Name0 = FALSE;

    ave_dmatrix(xyz, num, 3, ave_xyz);  /* default to all atoms */

    refmat = dmatrix(1, 4, 1, 3);
    identity_matrix(refmat, 3);  /* initialize it to identity matrix */
    if (args->pmi)  /* could update ave_xyz[] using only selected atoms */
        pdb_pmi(num, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, args->pmi, ave_xyz,
                refmat, rotmat);

    if (args->to_center)
        cpxyz(ave_xyz, rotmat[4]);
    move_position(xyz, num, 3, rotmat[4]);

    if (args->pmi) {  /* three views */
        tmprot = dmatrix(1, 3, 1, 3);
        multi_matrix(xyz, num, 3, rotmat, 3, 3, nxyz);  /* default: most extended */
        write_pdb(num, AtomName, ResName, ChainID, ResSeq, nxyz, Miscs, VIEW1_FILE);

        pmi_view(rotmat, tmprot, 3, 1, 2);  /* top view */
        multi_matrix(xyz, num, 3, tmprot, 3, 3, nxyz);
        write_pdb(num, AtomName, ResName, ChainID, ResSeq, nxyz, Miscs, VIEW2_FILE);

        pmi_view(rotmat, tmprot, 2, 3, 1);  /* landscape view */
        multi_matrix(xyz, num, 3, tmprot, 3, 3, nxyz);
        write_pdb(num, AtomName, ResName, ChainID, ResSeq, nxyz, Miscs, VIEW3_FILE);
        free_dmatrix(tmprot, 1, 3, 1, 3);
    }

    multi_matrix(xyz, num, 3, rotmat, 3, 3, nxyz);
    cmb_transformation(refmat, rotmat);
    free_dmatrix(refmat, 1, 4, 1, 3);

    write_pdb(num, AtomName, ResName, ChainID, ResSeq, nxyz, Miscs, args->outfile);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);

    free_cmatrix(Gvars.AtomName0, 1, DUMMY, 0, DUMMY);
    free_cmatrix(Gvars.ResName0, 1, DUMMY, 0, DUMMY);

    free_dmatrix(nxyz, 1, num, 1, 3);

}

static void rot_alcfile(double **rotmat, struct_args * args)
{
    char **AtomName;
    double **nxyz, **xyz, **refmat, **tmprot;
    long nbond, num, *ibase, **linkage;

    get_alc_nums(args->inpfile, &num, &nbond);
    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 1, 2);
    nxyz = dmatrix(1, num, 1, 3);

    read_alc(args->inpfile, &num, &nbond, AtomName, xyz, ibase, linkage);

    if (args->to_center)
        ave_dmatrix(xyz, num, 3, rotmat[4]);
    move_position(xyz, num, 3, rotmat[4]);

    refmat = dmatrix(1, 4, 1, 3);
    identity_matrix(refmat, 3);  /* initialize it to identity matrix */

    tmprot = dmatrix(1, 3, 1, 3);
    if (args->pmi) {  /* use all xyz coordinates */
        reorient_xyz(num, xyz, NULL, NULL, refmat);

        pmi_rotmat(num, xyz, rotmat);  /* default: most extended */
        multi_matrix(xyz, num, 3, rotmat, 3, 3, nxyz);
        write_alc(num, nbond, AtomName, nxyz, ibase, linkage, VIEW1_FILE);

        pmi_view(rotmat, tmprot, 3, 1, 2);  /* top view: rotx(-90) * rotz(-90) */
        multi_matrix(xyz, num, 3, tmprot, 3, 3, nxyz);
        write_alc(num, nbond, AtomName, nxyz, ibase, linkage, VIEW3_FILE);

        pmi_view(rotmat, tmprot, 2, 3, 1);  /* landscape view */
        multi_matrix(xyz, num, 3, tmprot, 3, 3, nxyz);
        write_alc(num, nbond, AtomName, nxyz, ibase, linkage, VIEW2_FILE);
    }
    free_dmatrix(tmprot, 1, 3, 1, 3);

    cmb_transformation(refmat, rotmat);
    free_dmatrix(refmat, 1, 4, 1, 3);

    multi_matrix(xyz, num, 3, rotmat, 3, 3, nxyz);

    write_alc(num, nbond, AtomName, nxyz, ibase, linkage, args->outfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_dmatrix(nxyz, 1, num, 1, 3);
}

/* rotate a structure [PDB or ALCHEMY] to adjust its orientation
   ---------------------------------------------------------------
   rotating or translating a PDB file in RasMol and then write it
   in PDB format does NOT change its coordinates. One way around
   is to "write molscript MOL-FILE", and then extract the rotation
   information from it, which looks like this:

   by rotation x 180.0
   by rotation z 105.669
   by rotation y -9.26690
   by rotation x 92.6453;

   the rotation matrix used here will be (note the negative values):
   rotx(-180)*rotz(-105.669)......

   to have a specific control over the orientation of a structure,
   either in PDB or ALCHEMY format, user may also need to create
   such a file by hand (with "by rotation" prefix).

   fnum = 0 means neither transformation matrix or rotation file
        = 1 means rotation file as from RasMol
        = 2 means transformation matrix
   --------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    double **rotmat;
    long num;
    struct_args args;

    set_my_globals(argv[0]);

    rotate_mol_cmdline(argc, argv, &args);

    if (args.fnum) {
        args.to_center = FALSE;
        args.pmi = 0;
        if (strcmp(args.rotfile, ROTMAT_FILE))
            remove_file(ROTMAT_FILE);
    } else
        remove_file(ROTMAT_FILE);
    remove_file(VIEW1_FILE);
    remove_file(VIEW2_FILE);
    remove_file(VIEW3_FILE);

    rotmat = dmatrix(1, 4, 1, 3);
    identity_matrix(rotmat, 3);
    if (args.fnum == 1)
        molscr_rotmat(args.rotfile, rotmat);
    else if (args.fnum == 2)
        read_rotmat(args.rotfile, rotmat);

    num = number_of_atoms(args.inpfile, 1, "*");
    (num) ? rot_pdbfile(num, rotmat, &args) : rot_alcfile(rotmat, &args);

    free_dmatrix(rotmat, 1, 4, 1, 3);

    clear_my_globals();

    return 0;
}
