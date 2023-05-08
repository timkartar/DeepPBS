#include "x3dna.h"

#define SIMPLE_BP_LONG_AXIS_RN9_YN1 2
#define SIMPLE_STEP_HELICAL_PARS 4

typedef struct {
    char torsion[BUF512];
    long istart;
    long istep;
    long icnt;
    long waters;
    long bz;  /* for B-Z junction */
    long ring;  /* for base ring center and normal vector */
    long simple_pars;  /* simplified base-pair (YC6--RC8) and step (C1'--C1') parameters */
    long abi;  /* Stephen Harvey's A/B index */
} struct_args;

static void get_ring_center_normal(long ds, long num_bp, long **pair_num, char **bp_seq,
                                   long **seidx, long *RY, char **AtomName,
                                   char **ResName, char *ChainID, long *ResSeq,
                                   char **Miscs, double **xyz, FILE * fp)
{
    static char *RingAtom[] = { RA_LIST };
    char BDIR[BUF512], idmsg[BUF512], sidmsg[BUF512], spdb[BUF512];
    char *sChainID, **sAtomName, **sResName, **sMiscs, *fmt = " %9.3f";

    double orgi[4], cxyz[4], nvec[4], rmsd;
    double **eRing_xyz, **fitted_xyz, **sRing_xyz, **sxyz, **R;

    long i, ib, ie, j, k, rnum, RingAtom_num, RA_NUM;
    long exp_katom, nmatch, snum, std_katom;
    long *sResSeq;

    get_BDIR(BDIR, "Atomic_A.pdb");
    RA_NUM = sizeof RingAtom / sizeof RingAtom[0];

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
        print_sep(fp, '*', 76);
        if (ds == 2)
            fprintf(fp, "Strand %s base ring center (Ox, Oy, Oz) and normal vector"
                    " (Nx, Ny, Nz) in\n   the coordinate system"
                    " of the given structure\n\n", (i == 1) ? "I" : "II");
        else
            fprintf(fp, "Base ring center (Ox, Oy, Oz) and normal vector (Nx, Ny, Nz) in"
                    "\n   the coordinate system of the given structure\n\n");

        fprintf(fp, "    base        Ox        Oy        Oz        Nx        Ny        Nz\n");

        for (j = 1; j <= num_bp; j++) {
            rnum = pair_num[i][j];
            ib = seidx[rnum][1];
            ie = seidx[rnum][2];
            get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);

            RingAtom_num = (RY[rnum] == 1) ? RA_NUM : RA_NUM - 3;
            set_std_base_pdb(BDIR, FALSE, bp_seq[i][j], spdb);
            snum = read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz,
                            sMiscs, 1, "*");
            sprintf(sidmsg, "in standard base: %s", spdb);

            nmatch = 0;
            for (k = 0; k < RingAtom_num; k++) {
                exp_katom = find_1st_atom(RingAtom[k], AtomName, ib, ie, idmsg);
                std_katom = find_1st_atom(RingAtom[k], sAtomName, 1, snum, sidmsg);
                if (exp_katom && std_katom) {
                    ++nmatch;
                    cpxyz(xyz[exp_katom], eRing_xyz[nmatch]);
                    cpxyz(sxyz[std_katom], sRing_xyz[nmatch]);
                }
            }

            rmsd = ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, orgi);
            UNUSED_PARAMETER(rmsd);

            ave_dmatrix(eRing_xyz, nmatch, 3, cxyz);
            for (k = 1; k <= 3; k++)  /* column-wise */
                nvec[k] = R[k][3];

            fprintf(fp, " %4ld %c   ", j, bp_seq[i][j]);
            for (k = 1; k <= 3; k++)
                fprintf(fp, fmt, cxyz[k]);  /* center */
            for (k = 1; k <= 3; k++)
                fprintf(fp, fmt, nvec[k]);  /* base normal */
            fprintf(fp, "\n");
        }
    }

    free_dmatrix(eRing_xyz, 1, RA_NUM, 1, 3);
    free_dmatrix(sRing_xyz, 1, RA_NUM, 1, 3);
    free_cmatrix(sAtomName, 1, NUM_RESIDUE_ATOMS, 0, 4);
    free_cmatrix(sResName, 1, NUM_RESIDUE_ATOMS, 0, 3);
    free_cvector(sChainID, 1, NUM_RESIDUE_ATOMS);
    free_lvector(sResSeq, 1, NUM_RESIDUE_ATOMS);
    free_dmatrix(sxyz, 1, NUM_RESIDUE_ATOMS, 1, 3);
    free_cmatrix(sMiscs, 1, NUM_RESIDUE_ATOMS, 0, NMISC);
    free_dmatrix(fitted_xyz, 1, RA_NUM, 1, 3);
    free_dmatrix(R, 1, 3, 1, 3);
}

static void simple_bp_pars(double **orien, double **org, long NNy, long **c6_c8,
                           long **chi, long num_bp, char **bp_seq, long *WC_info,
                           double **xyz, FILE * fp)
{
    char wc, str[BUF512], tmp[BUF32];
    char *bstr = "      ----", *fmt = "%10.2f";
    long i, j, k, Ly, Ry, num = 0;
    double buckle, propeller, r, angle;
    double dx[4], dy[4], dz[4], pars[7];
    double x1[4], y1[4], z1[4], x2[4], y2[4], z2[4];
    double morg[4], o1[4], o2[4], dd[4];
    double **r1, **r2, **mst, **parcln;

    print_sep(fp, '-', 76);
    fprintf(fp, "Simple base-pair parameters based on %s vectors\n",
            NNy ? "RN9--YN1" : "RC8--YC6");
    fprintf(fp, "      bp        Shear    Stretch   Stagger    Buckle  Propeller  Opening");
    fprintf(fp, "  angle\n");

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);
    parcln = dmatrix(1, num_bp, 1, 6);

    for (i = 1; i <= num_bp; i++) {
        /* Here: r1 refers to the base on strand II; r2 for strand I base.
         * For NEGATIVE orientaion, the r1 y/z-axes are already flipped */
        refs_right_left(i, orien, org, r1, o1, r2, o2);
        bpstep_par(r1, o1, r2, o2, pars, mst, morg);

        mtx_2_x_y_z(r1, x1, y1, z1);
        mtx_2_x_y_z(r2, x2, y2, z2);
        mtx_2_x_y_z(mst, dx, dy, dz);  /* keep dz as is */

        angle = z1_z2_angle_in_0_to_90(z1, z2);
        wc = (WC_info[i] == 2) ? ' ' : '*';
        sprintf(str, "%c %4ld %c%c%c ", wc, i, bp_seq[1][i], bp_seq[0][i], bp_seq[2][i]);

        if (NNy) {  /* RN9--YN1 */
            k = (i - 1) * 4;
            Ly = chi[1][k + 3];
            Ry = chi[2][k + 3];
        } else {  /* RC8--YC6 */
            Ly = c6_c8[1][i];
            Ry = c6_c8[2][i];
        }

        if (Ly && Ry) {
            num++;

            ddxyz(xyz[Ry], xyz[Ly], dd);
            vec_orth(dd, dz);  /* corrected y-axis */
            cpxyz(dd, dy);
            cross(dy, dz, dx);

            ddxyz(o1, o2, dd);
            pars[1] = dot(dd, dx);
            pars[2] = dot(dd, dy);
            pars[3] = dot(dd, dz);

            buckle = vec_ang(z1, z2, dx);
            propeller = vec_ang(z1, z2, dy);

            r = angle / sqrt(buckle * buckle + propeller * propeller);
            pars[4] = r * buckle;
            pars[5] = r * propeller;

            pars[6] = vec_ang(y1, y2, dz);

            for (j = 1; j <= 6; j++) {
                sprintf(tmp, fmt, pars[j]);
                strcat(str, tmp);
                parcln[num][j] = pars[j];
            }

        } else {
            for (j = 1; j <= 6; j++)
                strcat(str, bstr);
        }

        fprintf(fp, "%s %6.1f\n", str, angle);
    }

    output_ave_std(num, parcln, 1, fmt, fp);

    free_dmatrix(r1, 1, DUMMY, 1, DUMMY);
    free_dmatrix(r2, 1, DUMMY, 1, DUMMY);
    free_dmatrix(mst, 1, DUMMY, 1, DUMMY);
    free_dmatrix(parcln, 1, DUMMY, 1, DUMMY);
}

static void simple_step_heli_pars(double **orien, double **org, long **chi, long num_bp,
                                  char **bp_seq, long *WC_info, long *bphlx, double **xyz,
                                  long hel_pars, FILE * fp)
{
    char wc, str[BUF512], tmp[BUF32];
    char *bstr = "      ----", *fmt = "%10.2f";
    long i, j, ip1, k, c1a, c1b, num = 0;
    double morg[4], o1[4], o2[4], dd[4], pars[7];
    double *p1, *p2, **r1, **r2, **mst, **Rotmat;
    double **parcln;

    print_sep(fp, '-', 76);
    if (hel_pars) {
        fprintf(fp, "Simple base-pair helical parameters based on consecutive C1'-C1' vectors\n");
        fprintf(fp, "      step       X-disp    Y-disp   h-Rise     Incl.       Tip   h-Twist\n");
    } else {
        fprintf(fp, "Simple base-pair step parameters based on consecutive C1'-C1' vectors\n");
        fprintf(fp, "      step       Shift     Slide      Rise      Tilt      Roll     Twist\n");
    }

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);
    Rotmat = dmatrix(1, num_bp, 0, 12);  /* [0] as an indicator */

    /* collect pair frames */
    for (i = 1; i <= num_bp; i++) {
        refs_right_left(i, orien, org, r1, o1, r2, o2);
        bpstep_par(r1, o1, r2, o2, pars, mst, morg);

        p1 = Rotmat[i];  /* as a shorthand */
        mst2orien(p1, 0, mst);  /* keep z-axis as is */
        cpxyz(morg, p1 + 9);  /* keep bp origin in columns 10-12 */

        k = (i - 1) * 4;
        c1a = chi[1][k + 2];
        c1b = chi[2][k + 2];
        if (c1a && c1b) {
            p1[0] = TRUE;  /* with the C1'--C1' vector */
            ddxyz(xyz[c1b], xyz[c1a], dd);  /* 2-->1, for y-axis */
            vec_orth(dd, p1 + 6);  /* corrected y-axis */
            cpxyz(dd, p1 + 3);
            cross(p1 + 3, p1 + 6, p1);  /* x-axis */
        }
    }

    parcln = dmatrix(1, num_bp, 1, 6);
    for (i = 1; i < num_bp; i++) {
        ip1 = i + 1;
        p1 = Rotmat[i];
        p2 = Rotmat[ip1];
        orien2mst(p1, 0, r1);
        orien2mst(p2, 0, r2);

        wc = (WC_info[i] == 2 && WC_info[ip1] == 2) ? ' ' : '*';
        sprintf(str, "%c %4ld %c%c/%c%c", wc, i, bp_seq[1][i], bp_seq[1][ip1],
                bp_seq[2][ip1], bp_seq[2][i]);

        if (p1[0] && p2[0] && !bphlx[i]) {  /* with two sets of C1'--C1' vectors */
            ddxyz(p1 + 9, p2 + 9, dd);
            if (dot(dd, p1 + 6) < 0)
                reverse_y_z_columns(r1);
            if (dot(dd, p2 + 6) < 0)
                reverse_y_z_columns(r2);

            num++;
            if (hel_pars)
                helical_par(r1, p1 + 9, r2, p2 + 9, pars, mst, morg);
            else
                bpstep_par(r1, p1 + 9, r2, p2 + 9, pars, mst, morg);
            for (j = 1; j <= 6; j++) {
                sprintf(tmp, fmt, pars[j]);
                strcat(str, tmp);
                parcln[num][j] = pars[j];
            }

        } else
            for (j = 1; j <= 6; j++)
                strcat(str, bstr);

        fprintf(fp, "%s\n", str);
    }

    output_ave_std(num, parcln, 2, fmt, fp);

    free_dmatrix(r1, 1, DUMMY, 1, DUMMY);
    free_dmatrix(r2, 1, DUMMY, 1, DUMMY);
    free_dmatrix(mst, 1, DUMMY, 1, DUMMY);
    free_dmatrix(Rotmat, 1, DUMMY, 0, DUMMY);
    free_dmatrix(parcln, 1, DUMMY, 1, DUMMY);
}

static void check_simple_parameters(long simple, double **orien, double **org,
                                    long **c6_c8, long **chi, long num_bp, char **bp_seq,
                                    long *WC_info, long *bphlx, double **xyz, FILE * fp)
{
    long i, num, NNy, hel_pars, k = 0;

    if (!simple)
        return;

    for (i = 1; i <= num_bp; i++)
        if (WC_info[i] == 2)
            k++;
    num = num_bp - k;  /* number of non-WC base-pairs */

    print_sep(fp, '*', 76);
#include "simple-pars.inc"

    fprintf(fp, "This structure contains %ld non-Watson-Crick (with leading *) base pair(s)\n",
            num);

    NNy = simple & SIMPLE_BP_LONG_AXIS_RN9_YN1;  /* RN9--YN1 as (long) y-axis */
    simple_bp_pars(orien, org, NNy, c6_c8, chi, num_bp, bp_seq, WC_info, xyz, fp);

    simple_step_heli_pars(orien, org, chi, num_bp, bp_seq, WC_info, bphlx, xyz, FALSE, fp);

    hel_pars = simple & SIMPLE_STEP_HELICAL_PARS;  /* also derive/output helical parameters */
    if (hel_pars)
        simple_step_heli_pars(orien, org, chi, num_bp, bp_seq, WC_info, bphlx, xyz, TRUE, fp);
}

static void process_str(char *inpfile, struct_args * args)
{
    char pdbfile[BUF512], outfile[BUF512], **nt_info;
    char *ChainID, *bseq, **AtomName, **bp_seq, **ResName, **Miscs;
    double twist_p, twist_n, *mst_org, *mst_orgH, *mst_orien, *mst_orienH;
    double **org, **orien, **twist_rise, **xyz, **nt_bb_torsion;
    long i, ip, bbexist = 0, ds, hetatm = 0, parallel = 0, nbpm1, num;
    long num_bp, num_residue, str_type = 0;  /* for duplex only */
    long *idx, *WC_info, *ResSeq, *bphlx, *RY;
    long **c6_c8, **chi, **pair_num, **phos;
    long **seidx, **sugar, **o3p_brk, **htm_water;
    FILE *fp;

    /* get PDB file and pairing information */
    pair_num = read_input(inpfile, pdbfile, outfile, &ds, &num_bp, &ip, &hetatm);

    bphlx = lvector(1, num_bp);  /* helix break marker */
    for (i = 1; i <= num_bp; i++)
        bphlx[i] = args->icnt ? 0 : pair_num[ds + 1][i];

    fp = open_file(outfile, "w");

    /* read in PDB file */
    num = number_of_atoms(pdbfile, hetatm, Gvars.misc_pars.alt_list);
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(pdbfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs,
             hetatm, Gvars.misc_pars.alt_list);

    idx = lvector(1, num);
    atom_idx(num, AtomName, NULL, idx);

    /* get the numbering information of each residue */
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

    /* base-pairing residue number checking */
    pair_checking(ip, ds, num_residue, pdbfile, &num_bp, pair_num);

    /* check strand direction: O3'-P linkage, parallel etc */
    o3p_brk = lmatrix(1, ds, 1, num_bp);
    drct_checking(ds, num_bp, pair_num, seidx, AtomName, xyz, &parallel, &bbexist, o3p_brk, fp);

    /* get base or base-pair sequence */
    bp_seq = cmatrix(0, ds, 1, num_bp);  /* need to be kept */
    RY = lvector(1, num_residue);  /* simplified */
    get_bpseq(ds, num_bp, pair_num, seidx, AtomName, ResName, ChainID,
              ResSeq, Miscs, xyz, bp_seq, RY);

    /* for base overlap */
    bseq = cvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);  /* RY re-calculated */

    print_header(ds, num_bp, num, pdbfile, fp);

    /* get atom list for each residue */
    phos = lmatrix(1, ds + 4, 1, num_bp);  /* also include O1P & O2P atoms */
    c6_c8 = lmatrix(1, ds, 1, num_bp);
    sugar = lmatrix(1, ds, 1, num_bp * 5);
    chi = lmatrix(1, ds, 1, num_bp * 4);
    atom_list(ds, num_bp, pair_num, seidx, RY, bp_seq, AtomName, ResName,
              ChainID, ResSeq, Miscs, phos, c6_c8, sugar, chi);

    nt_info = cmatrix(1, num_residue, 0, BUF32);
    populate_nt_info(num_residue, seidx, ResName, ChainID, ResSeq, Miscs, bseq, nt_info);

    /* get the reference frame for each base */
    org = dmatrix(1, ds, 1, num_bp * 3);
    orien = dmatrix(1, ds, 1, num_bp * 9);
    WC_info = lvector(1, num_bp);
    ref_frames(ds, num_bp, pair_num, bp_seq, seidx, RY, AtomName, ResName, ChainID,
               ResSeq, Miscs, xyz, fp, orien, org, WC_info, &str_type, 0, o3p_brk);

    nt_bb_torsion = dmatrix(1, num_residue, 1, 6);
    get_nt_bb_torsion(nt_bb_torsion, num_residue, seidx, RY, AtomName, ResName,
                      ChainID, ResSeq, Miscs, xyz);

    if (ds == 2)  /* H-Bond information */
        hb_information(num_bp, pair_num, bp_seq, seidx, idx, AtomName, xyz, WC_info, fp);

    base_overlap(ds, num_bp, num, num_residue, pair_num, RY, bp_seq, seidx, AtomName, xyz, idx, orien, org, fp);  /* base overlap area in A^2 */

    /* get and print out parameters */
    nbpm1 = num_bp - 1;
    mst_orien = dvector(1, nbpm1 * 9);
    mst_org = dvector(1, nbpm1 * 3);
    mst_orienH = dvector(1, nbpm1 * 9);
    mst_orgH = dvector(1, nbpm1 * 3);
    twist_rise = dmatrix(1, nbpm1, 1, 2);
    get_parameters(ds, num_bp, bp_seq, orien, org, WC_info, fp, twist_rise, mst_orien,
                   mst_org, mst_orienH, mst_orgH, bphlx, args->istart, args->istep,
                   args->bz, &str_type, pair_num, nt_info);

    htm_water = lmatrix(1, 4, 0, num);  /* HETATM and water index */
    init_htm_water(args->waters, num, num_residue, idx, htm_water);
    identify_htw(num_residue, seidx, RY, AtomName, ResName, ChainID, ResSeq,
                 Miscs, xyz, htm_water);

    write_mst(ds, num_bp, pair_num, bp_seq, mst_orien, mst_org, seidx, AtomName,
              ResName, ChainID, ResSeq, xyz, Miscs, htm_water, twist_rise, STACK_FILE);
    write_mst(ds, num_bp, pair_num, bp_seq, mst_orienH, mst_orgH, seidx, AtomName,
              ResName, ChainID, ResSeq, xyz, Miscs, htm_water, twist_rise, HSTACK_FILE);

    if (args->ring)
        get_ring_center_normal(ds, num_bp, pair_num, bp_seq, seidx, RY, AtomName, ResName,
                               ChainID, ResSeq, Miscs, xyz, fp);

    if (ds == 2) {
        check_simple_parameters(args->simple_pars, orien, org, c6_c8, chi, num_bp, bp_seq,
                                WC_info, bphlx, xyz, fp);

        get_mtwist(nbpm1, bphlx, WC_info, twist_rise, &twist_p, &twist_n);

        str_classify(twist_p, twist_n, str_type, parallel, num_bp, fp);
        lambda_d3(num_bp, bp_seq, chi, c6_c8, xyz, fp);

        if (bbexist) {
            print_PP(parallel, twist_rise, num_bp, bp_seq, phos, mst_orien, mst_org,
                     mst_orienH, mst_orgH, xyz, WC_info, bphlx, args->abi, chi, fp);
            groove_width(parallel, num_bp, bp_seq, phos, xyz, bphlx, fp);
        }
        other_pars(num_bp, bp_seq, bphlx, orien, org);
    }

    /* global analysis */
    global_analysis(ds, num_bp, num, bp_seq, chi, phos, xyz, fp);

    if (bbexist) {
        backbone_torsion(ds, num_bp, pair_num, bp_seq, sugar, chi, xyz, nt_bb_torsion, fp);
        p_c1_dist(ds, num_bp, bp_seq, phos, chi, xyz, bphlx, fp);
        helix_radius(ds, num_bp, bp_seq, orien, org, phos, chi, xyz, bphlx, fp);
    }
    get_helix_axis(ds, num_bp, bp_seq, orien, org, bphlx, fp);

    /* free allocated vectors & matrices */
    free_lmatrix(pair_num, 1, ds + 1, 1, num_bp);
    free_lvector(bphlx, 1, num_bp);
    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_lmatrix(o3p_brk, 1, ds, 1, num_bp);
    free_cmatrix(bp_seq, 0, ds, 1, num_bp);
    free_lvector(RY, 1, num_residue);
    free_cvector(bseq, 1, num_residue);
    free_lmatrix(phos, 1, ds + 4, 1, num_bp);
    free_lmatrix(c6_c8, 1, ds, 1, num_bp);
    free_dmatrix(nt_bb_torsion, 1, num_residue, 1, 6);
    free_lmatrix(sugar, 1, ds, 1, num_bp * 5);
    free_lmatrix(chi, 1, ds, 1, num_bp * 4);
    free_dmatrix(orien, 1, ds, 1, num_bp * 9);
    free_dmatrix(org, 1, ds, 1, num_bp * 3);
    free_lvector(WC_info, 1, num_bp);
    free_dvector(mst_orien, 1, nbpm1 * 9);
    free_dvector(mst_org, 1, nbpm1 * 3);
    free_dvector(mst_orienH, 1, nbpm1 * 9);
    free_dvector(mst_orgH, 1, nbpm1 * 3);
    free_dmatrix(twist_rise, 1, nbpm1, 1, 2);
    free_lmatrix(htm_water, 1, 4, 0, num);
    free_lvector(idx, 1, num);
    free_cmatrix(nt_info, 1, num_residue, 0, BUF32);

    close_file(fp);
}

static void calculate_torsions(char *outfile, char *pdbfile)
{
    char BDIR[BUF512], **nt_info;
    char *ChainID, *bseq, **AtomName, **ResName, **Miscs;
    double **org, **orien, **xyz, **nt_torsion, **ss_Zp_Dp;
    long num, num_residue, hetatm = TRUE;
    long *ResSeq, *RY, **seidx, **nt_list;
    FILE *fp;

    fp = open_file(outfile, "w");

    num = number_of_atoms(pdbfile, hetatm, Gvars.misc_pars.alt_list);
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(pdbfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs,
             hetatm, Gvars.misc_pars.alt_list);

    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

    bseq = cvector(1, num_residue);
    RY = lvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);

    org = dmatrix(1, num_residue, 1, 3);
    orien = dmatrix(1, num_residue, 1, 9);
    get_BDIR(BDIR, "Atomic_A.pdb");
    base_frame(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID, ResSeq, Miscs,
               xyz, BDIR, orien, org);

    nt_info = cmatrix(1, num_residue, 0, BUF32);
    populate_nt_info(num_residue, seidx, ResName, ChainID, ResSeq, Miscs, bseq, nt_info);

    nt_list = lmatrix(1, num_residue, 0, BUF32);
    populate_nt_list(num_residue, seidx, RY, bseq, AtomName, xyz, nt_list);
    output_Borg_P_C1_C4(num_residue, org, xyz, nt_list, nt_info);

    nt_torsion = dmatrix(1, num_residue, 1, BUF32);
    get_nt_torsion(num_residue, org, xyz, nt_list, nt_torsion);

    ss_Zp_Dp = dmatrix(1, 2, 1, num_residue);
    get_ss_Zp_Dp(num_residue, org, orien, xyz, nt_list, ss_Zp_Dp);

    output_nt_torsion(num_residue, nt_info, nt_list, nt_torsion, ss_Zp_Dp, fp);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_lvector(RY, 1, num_residue);
    free_cvector(bseq, 1, num_residue);
    free_dmatrix(org, 1, num_residue, 1, 3);
    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(nt_torsion, 1, num_residue, 1, BUF32);
    free_lmatrix(nt_list, 1, num_residue, 0, BUF32);
    free_cmatrix(nt_info, 1, num_residue, 0, BUF32);
    free_dmatrix(ss_Zp_Dp, 1, 2, 1, num_residue);

    close_file(fp);
}

static void analyze_usage(void)
{
    help3dna_usage("analyze");
}

static void set_defaults(struct_args * args)
{
    strcpy(args->torsion, "");
    args->istart = 1;
    args->istep = 1;
    args->icnt = FALSE;
    args->waters = FALSE;
    args->bz = TRUE;  /* check for B-Z junction */
    args->ring = FALSE;
    args->simple_pars = TRUE;
    args->abi = FALSE;
}

static long derive_simple_pars_settings(char *option)
{
    long k = FALSE;

    if (lux_ncmatch(option, "no|false|off"))
        return k;

    k = TRUE;  /* default: RC8--YC6, no helical pars */

    if (lux_ncmatch(option, "n1|n9"))
        k |= SIMPLE_BP_LONG_AXIS_RN9_YN1;

    if (lux_ncmatch(option, "heli"))
        k |= SIMPLE_STEP_HELICAL_PARS;

    fprintf(stderr, "Setting for simple parameters: %ld\n", k);
    return k;
}

static void analyze_cmdline(int argc, char *argv[], struct_args * args)
{
    long i, j, k;

    set_defaults(args);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        if (lux_ncmatch(argv[i], "^--?t")) {
            get_strvalue(argv[i], args->torsion, FALSE);
            continue;
        }

        if (lux_ncmatch(argv[i], "^--?bz")) {
            args->bz = set_switch_default_true(argv[i]);
            continue;
        }

        if (lux_ncmatch(argv[i], "^--?ri")) {
            args->ring = set_switch_default_true(argv[i]);
            continue;
        }

        if (lux_ncmatch(argv[i], "^--?si")) {
            args->simple_pars = derive_simple_pars_settings(argv[i]);
            continue;
        }

        if (lux_ncmatch(argv[i], "^--?abi?")) {
            args->abi = set_switch_default_true(argv[i]);
            continue;
        }

        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'C')
                args->icnt = TRUE;
            else if (argv[i][j] == 'W')
                args->waters = TRUE;
            else if (argv[i][j] == 'S') {
                k = sscanf(argv[i], "-S=%ld,%ld", &args->istep, &args->istart);
                if (k == 2) {
                    if (args->istart < 0)
                        args->istart = -args->istart;  /* make it positive */
                } else if (k == 1)
                    args->istart = 1;  /* redundant: already initialized */
                else {
                    fprintf(stderr, "wrong format for setting step\n");
                    analyze_usage();
                }
                fprintf(stderr, "***start at %ld, with step size: %ld***\n", args->istart,
                        args->istep);
                break;  /* -s=istep,istart not combined with others */
            } else
                analyze_usage();
    }

    if (argc == i) {
        if (is_empty_string(args->torsion))
            process_str("stdin", args);
        else
            analyze_usage();
    } else {
        if (is_empty_string(args->torsion)) {  /* regular case */
            for (j = i; j < argc; j++) {
                if (strcmp(argv[j], "tmpfile"))  /* analyze multiple structures */
                    fprintf(stderr, "\n......Processing structure #%ld: <%s>......\n",
                            j - i + 1, argv[j]);
                process_str(argv[j], args);
            }
        } else {  /* torsion angle only */
            calculate_torsions(args->torsion, argv[i]);  /* only one PDB structure */
        }
    }
}

int main(int argc, char *argv[])
{
    struct_args args;
    time_t time0;

    time(&time0);

    /* clean up files with fixed name */
    remove_file(AUX_FILE);
    remove_file(BPSTEP_FILE);
    remove_file(HLXSTEP_FILE);
    remove_file(STACK_FILE);
    remove_file(HSTACK_FILE);
    remove_file(REF_FILE);
    remove_file(POC_FILE);
    remove_file(SEVEN_FILE);

    set_my_globals(argv[0]);
    analyze_cmdline(argc, argv, &args);
    clear_my_globals();

    print_used_time(time0);

    return 0;
}
