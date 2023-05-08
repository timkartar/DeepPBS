#include "x3dna.h"

static void anyhelix_usage(void)
{
    help3dna_usage("anyhelix");
}

/* get PDB data file name, output file name, and pairing information
   pair_num matrix is allocated here but de-allocated elsewhere
   The function is adapted from "read_input" of analyze, and is designed to
   handle any helix: single, double, triplex, etc, parallel or anti-parallel */
static long **read_pair(char *inpfile, char *pdbfile, char *outfile, long *num_bp,
                        long *hetatm, long *eqnum)
{
    char str[BUF512], strc[BUF512];
    char *item[MBASES];
    long i, inum, j, k, num_bases = 0, nitem;
    long **pair_num;
    FILE *fp;

    fp = open_file(inpfile, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%s", pdbfile) != 1)
        fatal("can not read pdbfile name: %s\n", pdbfile);
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%s", outfile) != 1)
        fatal("can not read outfile name: %s\n", outfile);
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", &num_bases) != 1)
        fatal("can not read number of bases per-pair-set\n");
    if (num_bases <= 0)
        fatal("illegal number (%ld <= 0) of bases per-pair-set\n", num_bases);
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", num_bp) != 1)
        fatal("can not read number of base-pairs\n");
    if (*num_bp <= 0)
        fatal("illegal number (%ld <= 0) of base-pairs\n", *num_bp);
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld %ld", &k, hetatm) != 2)
        fatal("can not read integer pair numbering/hetero indicator\n");

    *eqnum = 1;
    pair_num = lmatrix(1, *num_bp, 0, MBASES);
    for (i = 1; i <= *num_bp; i++) {
        if (fgets(str, sizeof str, fp) == NULL)
            fatal("can not read pair list\n");
        strcpy(strc, str);  /* itemize will damage str */
        nitem = itemize(str, item, MBASES) + 1;
        inum = 0;
        for (j = 0; j < nitem; j++) {
            if (inum >= num_bases || sscanf(item[j], "%ld", &k) != 1)
                break;
            if (k <= 0)
                fatal("residue number out of range (%ld <= 0)!\n", k);
            pair_num[i][++inum] = k;
        }
        if (!j)
            fatal("containing line with no base residue numbers\n");
        if (inum != num_bases) {
            *eqnum = 0;
            fprintf(stderr, "%s           ==> has %ld instead of %ld bases\n",
                    strc, inum, num_bases);
        }
        pair_num[i][0] = inum;
    }
    close_file(fp);

    return pair_num;
}

/* cross checking all possible pairs in a multiplet */
static void cross_pair(long num, long num_residue, long num_bp, long **pair_num,
                       miscPars * misc_pars, char *bseq, long **seidx, char **AtomName,
                       char **ResName, char *ChainID, long *ResSeq, char **Miscs,
                       double **xyz, double *chi, double **orien, double **org,
                       double **NC1xyz, double *bp_orien, double *bp_org, long *RY,
                       long *idx, long **htm_water, FILE * fp)
{
    char wc[BUF512], str[BUF512], b1[BUF512], b2[BUF512], tmp[BUF512];
    double rtn_val[RTNNUM];
    long i, ia, ib, ik, inum, j, jr, k, kr;
    long bpid, ntot = 0, ivec[3], **ring_atom;
    FILE *mfp, *fp2;

    /* 1-9 ring atom index, 10 # of ring atoms, 11-19 first level */
    ring_atom = lmatrix(1, num_residue, 1, 19);
    ring_oidx(num, num_residue, RY, seidx, AtomName, xyz, idx, ring_atom);

    print_bp_crit(misc_pars, fp);
    fprintf(fp, "For each base-layer (a base-pair/triplet etc), detailed residue"
            " information\n    is given followed by a description of each possible"
            " base-pair in the layer,\n    its base-pair parameters (shear, stretch,"
            " stagger, buckle, propeller,\n    opening) and H-bonding information.\n");

    mfp = open_file(MUL_FILE, "w");
    fp2 = open_file(ALLP_FILE, "w");
    for (i = 1; i <= num_bp; i++) {  /* each possible pair */
        inum = pair_num[i][0];
        str[0] = '\0';
        for (j = 1; j <= inum; j++) {
            k = pair_num[i][j];
            kr = seidx[k][1];
            base_str(ChainID[kr], ResSeq[kr], Miscs[kr], ResName[kr], bseq[k], 1, b1);
            sprintf(tmp, "[%ld]%s%s", k, b1, (j == inum) ? "" : " + ");
            strcat(str, tmp);
        }
        fprintf(fp, "\n%4ld (%ld): %s\n", i, inum, str);
        fprintf(mfp, "%6s    %4ld\n", "MODEL ", i);
        fprintf(mfp, "REMARK    Section #%4.4ld %ld\n", i, inum);
        fprintf(mfp, "REMARK    %s\n", str);
        fprintf(mfp, "REMARK    %s\n", Gvars.X3DNA_VER);
        pair2mst(inum, pair_num[i], AtomName, ResName, ChainID, ResSeq, Miscs, xyz,
                 orien, org, seidx, bp_orien + (i - 1) * 9, bp_org + (i - 1) * 3,
                 htm_water, misc_pars, mfp);
        fprintf(mfp, "ENDMDL\n");
        ik = 0;
        for (j = 1; j < inum; j++) {
            for (k = j + 1; k <= inum; k++) {
                jr = pair_num[i][j];
                kr = pair_num[i][k];
                check_pair(jr, kr, bseq, seidx, xyz, NC1xyz, orien, org, idx, AtomName,
                           misc_pars, rtn_val, &bpid, ring_atom, 0);
                if (bpid) {  /* j & k paired */
                    bpid_wc_str(bpid, rtn_val[35], wc);
                    ia = seidx[jr][1];
                    ib = seidx[kr][1];
                    base_str(ChainID[ia], ResSeq[ia], Miscs[ia], ResName[ia], bseq[jr], 1, b1);
                    base_str(ChainID[ib], ResSeq[ib], Miscs[ib], ResName[ib], bseq[kr], 2, b2);
                    sprintf(tmp, "%s-%s-%s", b1, wc, b2);
                    fprintf(fp, "   %4ld (%ld): [%ld]%s[%ld]\n", ++ntot, ++ik, jr, tmp, kr);
                    print_pairinfo(jr, kr, bseq[jr], bseq[kr], rtn_val, chi,
                                   misc_pars, seidx, idx, AtomName, xyz, bseq, 0, fp);

                    fprintf(fp2, "%6s    %4ld\n", "MODEL ", ntot);
                    fprintf(fp2, "REMARK    Section #%4.4ld %s\n", ntot, tmp);
                    fprintf(fp2, "REMARK    %s\n", Gvars.X3DNA_VER);
                    ivec[1] = jr;
                    ivec[2] = kr;
                    pair2mst(2, ivec, AtomName, ResName, ChainID, ResSeq, Miscs, xyz,
                             orien, org, seidx, NULL, NULL, htm_water, misc_pars, fp2);
                    fprintf(fp2, "ENDMDL\n");
                }
            }
        }
    }

    free_lmatrix(ring_atom, 1, num_residue, 1, 19);
    close_file(mfp);
    close_file(fp2);
}

/* get the base set step in string */
static void stepstr(long i, long j, char *bseq, long **pair_num, char *str)
{
    long k, n;

    n = sprintf(str, "%4ld ", i);
    for (k = 1; k <= pair_num[i][0]; k++)
        str[n++] = bseq[pair_num[i][k]];
    if (j > 0) {
        str[n++] = '/';
        for (k = 1; k <= pair_num[j][0]; k++)
            str[n++] = bseq[pair_num[j][k]];
    }
    str[n] = '\0';
}

/* print local base-pair step and helical parameter for any helix */
static void print_ah_par(long ishel, long num_bp, char *bseq, long **pair_num,
                         double **param, FILE * fp)
{
    char *format = "%10.2f", lsp[BUF512], str[BUF512];
    double temp[7];
    long i, j, nbpm1, nlsp = 0;

    if (num_bp == 1)
        return;

    if (pair_num[1][0] > 3) {
        nlsp = (pair_num[1][0] - 3) * 2;
        for (i = 0; i < nlsp; i++)
            lsp[i] = ' ';
    }
    lsp[nlsp] = '\0';

    if (ishel)
        fprintf(fp, "%s      step       X-disp    Y-disp   h-Rise"
                "     Incl.       Tip   h-Twist\n", lsp);
    else
        fprintf(fp, "%s      step       Shift     Slide      Rise"
                "      Tilt      Roll     Twist\n", lsp);
    nbpm1 = num_bp - 1;
    for (i = 1; i <= nbpm1; i++) {
        stepstr(i, i + 1, bseq, pair_num, str);
        fprintf(fp, "%s", str);
        for (j = 1; j <= 6; j++)
            fprintf(fp, format, param[i][j]);
        fprintf(fp, "\n");
    }

    if (nbpm1 >= 2) {
        fprintf(fp, "            ");
        print_sep(fp, '~', 60 + nlsp);
        fprintf(fp, "%s        ave.", lsp);
        ave_dmatrix(param, nbpm1, 6, temp);
        for (i = 1; i <= 6; i++)
            fprintf(fp, format, temp[i]);
        fprintf(fp, "\n");
        fprintf(fp, "%s        s.d.", lsp);
        std_dmatrix(param, nbpm1, 6, temp);
        for (i = 1; i <= 6; i++)
            fprintf(fp, format, temp[i]);
        fprintf(fp, "\n");
    }
}

/* write out a multiple base step with reference to its middle frame */
static void ah_write_mst(long i, long **pair_num, char *bseq, long **seidx, double **mst,
                         double *msto, char **AtomName, char **ResName, char *ChainID,
                         long *ResSeq, char **Miscs, double **xyz, long **htm_water,
                         miscPars * misc_pars, FILE * fp)
{
    char str[BUF512];
    double **xyz_residue;
    long ik, j, jr, k, tnum_res1, tnum_res2;
    long ivec0[BUF512], ivect[BUF512];
    long inum = 0;

    xyz_residue = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);

    ik = i + 1;
    stepstr(i, ik, bseq, pair_num, str);

    fprintf(fp, "%6s    %4ld\n", "MODEL ", i);
    fprintf(fp, "REMARK    Section #%4.4ld %s\n", i, str);
    fprintf(fp, "REMARK    %s\n", Gvars.X3DNA_VER);

    tnum_res1 = attached_residues(pair_num[i][0], pair_num[i], ivec0, seidx,
                                  xyz, htm_water, misc_pars);
    fprintf(fp, "REMARK    LOWER: %4ld", tnum_res1);
    for (j = 1; j <= tnum_res1; j++) {  /* lower base */
        ivect[j] = ivec0[j];
        fprintf(fp, "%4ld", j);
    }
    fprintf(fp, "\n");

    tnum_res2 = attached_residues(pair_num[ik][0], pair_num[ik], ivec0, seidx,
                                  xyz, htm_water, misc_pars);
    fprintf(fp, "REMARK    UPPER: %4ld", tnum_res2);
    for (j = 1; j <= tnum_res2; j++) {
        k = j + tnum_res1;
        ivect[k] = ivec0[j];  /* upper bases */
        fprintf(fp, "%4ld", k);
    }
    fprintf(fp, "\n");

    for (j = 1; j <= tnum_res1 + tnum_res2; j++) {
        jr = ivect[j];
        for (k = seidx[jr][1]; k <= seidx[jr][2]; k++)
            cpxyz(xyz[k], xyz_residue[k - seidx[jr][1] + 1]);
        change_xyz(0, msto, mst, seidx[jr][2] - seidx[jr][1] + 1, xyz_residue);
        pdb_record(seidx[jr][1], seidx[jr][2], &inum, 1, AtomName,
                   ResName, ChainID, ResSeq, xyz_residue, Miscs, fp);
    }
    fprintf(fp, "ENDMDL\n");

    free_dmatrix(xyz_residue, 1, NUM_RESIDUE_ATOMS, 1, 3);
}

/* get the local step and helical parameters of any helix */
static void ah_step_hel(long num_bp, long **pair_num, char *bseq, long **seidx,
                        double *bp_orien, double *bp_org, char **AtomName, char **ResName,
                        char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                        long **htm_water, FILE * fp)
{
    double hfoi[4], mfoi[4], o1[4], o2[4];
    double **bp_step_par, **bp_heli_par;
    double **mfi, **hfi, **r1, **r2;
    long i, nbpm1;
    FILE *fps, *fph;

    nbpm1 = num_bp - 1;
    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mfi = dmatrix(1, 3, 1, 3);
    hfi = dmatrix(1, 3, 1, 3);
    bp_step_par = dmatrix(1, nbpm1, 1, 6);
    bp_heli_par = dmatrix(1, nbpm1, 1, 6);

    fps = open_file(STACK_FILE, "w");
    fph = open_file(HSTACK_FILE, "w");

    for (i = 1; i <= nbpm1; i++) {
        refs_i_j(i, i + 1, bp_orien, bp_org, r1, o1, r2, o2);
        bpstep_par(r1, o1, r2, o2, bp_step_par[i], mfi, mfoi);
        helical_par(r1, o1, r2, o2, bp_heli_par[i], hfi, hfoi);

        ah_write_mst(i, pair_num, bseq, seidx, mfi, mfoi, AtomName, ResName, ChainID,
                     ResSeq, Miscs, xyz, htm_water, &Gvars.misc_pars, fps);
        ah_write_mst(i, pair_num, bseq, seidx, hfi, hfoi, AtomName, ResName, ChainID,
                     ResSeq, Miscs, xyz, htm_water, &Gvars.misc_pars, fph);
    }

    print_sep(fp, '*', 76);
    fprintf(fp, "Local base-pair step parameters\n");
    print_ah_par(0, num_bp, bseq, pair_num, bp_step_par, fp);

    print_sep(fp, '*', 76);
    fprintf(fp, "Local base-pair helical parameters\n");
    print_ah_par(1, num_bp, bseq, pair_num, bp_heli_par, fp);

    close_file(fps);
    close_file(fph);

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mfi, 1, 3, 1, 3);
    free_dmatrix(hfi, 1, 3, 1, 3);
    free_dmatrix(bp_step_par, 1, nbpm1, 1, 6);
    free_dmatrix(bp_heli_par, 1, nbpm1, 1, 6);
}

/* get and print P, O4' and C1' atom radius from global helix in Raster3D format */
static void get_ah_poc_radius(long num_bp, long **pair_num, long **seidx, char **AtomName,
                              double **xyz, double *haxis, double *hstart, double *hend,
                              double *rave, double *rstd)
{
    static char *poc_atom[3] = { " P  ", " O4'", " C1'" };
    long i, ib, ie, ik, j, jr, k, num_ple, poc_num[3];
    double **poc_radius;

    for (i = 0; i <= 2; i++) {
        poc_num[i] = 0;
        rave[i] = rstd[i] = 0.0;
    }

    num_ple = pair_num[1][0] * num_bp;  /* maximum possible number */
    poc_radius = dmatrix(0, 2, 1, num_ple);
    for (i = 1; i <= num_bp; i++)
        for (j = 1; j <= pair_num[i][0]; j++) {
            jr = pair_num[i][j];
            ib = seidx[jr][1];
            ie = seidx[jr][2];
            for (k = 0; k <= 2; k++) {
                ik = find_1st_atom(poc_atom[k], AtomName, ib, ie, "");
                if (ik)
                    poc_radius[k][++poc_num[k]] = a_hlxdist(ik, xyz, haxis, hstart);
            }
        }

    for (i = 0; i <= 2; i++)
        if (poc_num[i]) {
            rave[i] = ave_dvector(poc_radius[i], poc_num[i]);
            rstd[i] = std_dvector(poc_radius[i], poc_num[i]);
        }
    print_poc_r3d(rave, hstart, hend);

    free_dmatrix(poc_radius, 0, 2, 1, num_ple);
}

/* get twist and rise based on C1'-C1' vectors: should check if C1' exists */
static void c1_twist_rise(long i, long nb, double *mz, double **NC1xyz, long **pair_num, FILE * fp)
{
    long i1, i2, j, j1_osx, j2;
    double d1[4], d2[4], *ctwist, *crise;

    ctwist = dvector(1, nb);
    crise = dvector(1, nb);

    vec_norm(mz);
    fprintf(fp, "     ");  /* rise */
    for (j = 1; j <= nb; j++) {
        i1 = pair_num[i][j];
        i2 = pair_num[i + 1][j];
        ddxyz(NC1xyz[i1] + 3, NC1xyz[i2] + 3, d1);
        crise[j] = dot(d1, mz);
        fprintf(fp, "%7.2f", crise[j]);
    }
    d1[1] = ave_dvector(crise, nb);
    d1[2] = (nb > 1) ? std_dvector(crise, nb) : 0.0;
    fprintf(fp, "  ==> %7.2f(%.2f) (Angstroms)\n", d1[1], d1[2]);

    fprintf(fp, "     ");  /* twist */
    for (j = 1; j <= nb; j++) {
        i1 = pair_num[i][j];
        i2 = pair_num[i + 1][j];
        if (j == nb) {  /* round to beginning */
            j1_osx = pair_num[i][1];
            j2 = pair_num[i + 1][1];
        } else {
            j1_osx = pair_num[i][j + 1];
            j2 = pair_num[i + 1][j + 1];
        }
        ddxyz(NC1xyz[i1] + 3, NC1xyz[j1_osx] + 3, d1);
        ddxyz(NC1xyz[i2] + 3, NC1xyz[j2] + 3, d2);
        if (nb > 1 || j < nb) {
            ctwist[j] = vec_ang(d1, d2, mz);
            fprintf(fp, "%7.1f", ctwist[j]);
        }
    }
    d1[1] = (nb > 1) ? ave_dvector(ctwist, nb) : 0.0;
    d1[2] = (nb > 1) ? std_dvector(ctwist, nb) : 0.0;
    fprintf(fp, "  %s==> %7.1f(%.1f) (degrees)\n", (nb < 2) ? "       " : "", d1[1], d1[2]);

    free_dvector(ctwist, 1, nb);
    free_dvector(crise, 1, nb);
}

/* get the global helical axis of anyhelix */
static void get_ah_helix(long num_bp, char *bseq, long **pair_num, long num, long **idxCN,
                         long **seidx, char **AtomName, double **xyz, double **NC1xyz,
                         double *bp_org, FILE * fp)
{
    char str[BUF512];
    double hrise, std_rise, rave[4], rstd[4], haxis[4], hstart[4], hend[4];
    long i, ia, ib, j, k, nb, ntot, nvec = 0;
    long C1b[MBASES], C1e[MBASES], **idx;

    nb = pair_num[1][0];
    ntot = nb * 2 * (num_bp - 1);
    idx = lmatrix(1, ntot, 1, 2);  /* beginning & end index */

    for (i = 1; i <= num_bp - 1; i++)
        for (j = 1; j <= nb; j++)
            for (k = 1; k <= 2; k++) {
                ia = idxCN[pair_num[i][j]][k];  /* C/N # of i */
                ib = idxCN[pair_num[i + 1][j]][k];  /* C/N # of i + 1 */
                if (ia && ib) {  /* both exist */
                    idx[++nvec][1] = ia;
                    idx[nvec][2] = ib;
                }
            }
    for (i = 1; i <= nb; i++) {  /* assuming 1st & last C1* atoms exist */
        C1b[i] = idxCN[pair_num[1][i]][1];
        C1e[i] = idxCN[pair_num[num_bp][i]][1];
    }
    get_axis(nvec, idx, num, xyz, nb, C1b, C1e, &std_rise, &hrise, haxis, hstart, hend);
    free_lmatrix(idx, 1, ntot, 1, 2);
    if (hrise < EMPTY_CRITERION)
        return;  /* no helix defined */

    print_sep(fp, '*', 76);
    fprintf(fp, "Global linear helical axis defined by equivalent C1'"
            " and RN9/YN1 atom pairs\n");
    fprintf(fp, "Deviation from regular linear helix: %.2f(%.2f)\n", hrise, std_rise);

    if (std_rise > Gvars.misc_pars.std_curved)
        return;

    get_ah_poc_radius(num_bp, pair_num, seidx, AtomName, xyz, haxis, hstart, hend, rave, rstd);
    fprintf(fp, "Helix:  %8.3f%8.3f%8.3f\n", haxis[1], haxis[2], haxis[3]);
    fprintf(fp, "HETATM 9998  XS    X X 999    %8.3f%8.3f%8.3f\n",
            hstart[1], hstart[2], hstart[3]);
    fprintf(fp, "HETATM 9999  XE    X X 999    %8.3f%8.3f%8.3f\n", hend[1], hend[2], hend[3]);
    fprintf(fp, "Average and standard deviation of helix radius:\n");
    fprintf(fp, "      P: %.2f(%.2f), O4': %.2f(%.2f),  C1': %.2f(%.2f)\n",
            rave[0], rstd[0], rave[1], rstd[1], rave[2], rstd[2]);

    fprintf(fp, "\nGlobal rise and twist based on C1'-C1' vectors\n"
            "Number in [] is the projection of the vector connecting neighboring mean\n"
            "   base-layer centers onto the global helix defined above\n"
            "Numbers following ==> are mean value and standard deviation\n");
    for (i = 1; i < num_bp; i++) {
        stepstr(i, i + 1, bseq, pair_num, str);
        ddxyz(bp_org + (i - 1) * 3, bp_org + i * 3, rave);
        fprintf(fp, "%s [%5.2f Angstroms]\n", str, dot(haxis, rave));
        c1_twist_rise(i, nb, haxis, NC1xyz, pair_num, fp);
    }
}

/* print out the reference frame of each base set for "frame_mol" */
static void print_ah_frame(long num_bp, long **pair_num, char *bseq, double *bp_orien,
                           double *bp_org)
{
    char *format = "%10.4f%10.4f%10.4f\n";
    long i, ik, ioffset3, ioffset9, j;
    FILE *fp;

    fp = open_file(REF_FILE, "w");

    fprintf(fp, "%5ld bases\n", num_bp);
    for (i = 1; i <= num_bp; i++) {
        fprintf(fp, "... %5ld ", i);
        for (j = 1; j <= pair_num[i][0]; j++)
            fprintf(fp, "%c", bseq[pair_num[i][j]]);
        fprintf(fp, " ...\n");
        ioffset3 = (i - 1) * 3;
        fprintf(fp, format, bp_org[ioffset3 + 1], bp_org[ioffset3 + 2], bp_org[ioffset3 + 3]);
        ioffset9 = (i - 1) * 9;
        for (j = 1; j <= 3; j++) {
            ik = (j - 1) * 3;
            fprintf(fp, format, bp_orien[ioffset9 + ik + 1],
                    bp_orien[ioffset9 + ik + 2], bp_orien[ioffset9 + ik + 3]);
        }
    }
    close_file(fp);
}

/* reverse x- and z-axes for a layer of bases (as in handling Z-DNA) */
static void reverse_xz(long idx, double **orien, long **pair_num, long *s1dir)
{
    long i, k;

    s1dir[idx] = 1;
    for (i = 1; i <= pair_num[idx][0]; i++) {
        k = pair_num[idx][i];
        negate_xyz(orien[k]);  /* x-axis */
        negate_xyz(orien[k] + 6);  /* z-axis */
    }
}

/* to account for flipped over bases which make z-axes reverse their directions (as in
   ud0013). Strand 1's z-axis should always be along its 5'--->3' direction: otherwise
   rotate about y-axis by 180 degrees: as in Z-DNA */
static void stnd1_up(double **orien, double **org, double *chi, long num_bp,
                     long **pair_num, FILE * fp)
{
    long i, j, *s1dir;
    double dxyz[4];

    if (num_bp == 1)
        return;

    s1dir = lvector(1, num_bp);

    /* check the first base on strand I */
    ddxyz(org[pair_num[1][1]], org[pair_num[2][1]], dxyz);
    if (dot(dxyz, &orien[pair_num[1][1]][6]) < 0.0)
        reverse_xz(1, orien, pair_num, s1dir);

    /* all following strand 1 bases follow base1 */
    for (i = 2; i <= num_bp; i++) {
        ddxyz(org[pair_num[i - 1][1]], org[pair_num[i][1]], dxyz);
        if (dot(&orien[pair_num[i][1]][6], &orien[pair_num[i - 1][1]][6]) < 0.0 &&
            dval_in_range(dot(dxyz, &orien[pair_num[i - 1][1]][6]), 1.5, 7.5))
            reverse_xz(i, orien, pair_num, s1dir);
    }

    for (i = 1; i <= num_bp; i++) {
        j = pair_num[i][1];
        fprintf(fp, "Residue %ld has chi = %.2f: %s\n", j, chi[j],
                (s1dir[i]) ? "x- and z-axes REVERSED" : "");
    }

    free_lvector(s1dir, 1, num_bp);
}

/* get bend angle, local twist and rise based on C1-C1 vectors:
 * equal bases per layer pair_num[i][0] = pair_num[1][0] */
static void local_c1par(long num_bp, char *bseq, double *bp_orien, double **NC1xyz,
                        long **pair_num, FILE * fp)
{
    char str[BUF512];
    double z[4];
    long i;

    print_sep(fp, '*', 76);
    fprintf(fp, "Local bend angle in [], rise and twist based on C1'-C1' vectors\n"
            "Numbers following ==> are mean value and standard deviation\n");
    for (i = 1; i < num_bp; i++) {
        stepstr(i, i + 1, bseq, pair_num, str);
        fprintf(fp, "%s [%5.1f degrees]\n", str,
                magang(bp_orien + i * 9 - 3, bp_orien + i * 9 + 6));
        sumxyz(bp_orien + i * 9 + 6, bp_orien + (i - 1) * 9 + 6, z);  /* mean z-axis */
        c1_twist_rise(i, pair_num[i][0], z, NC1xyz, pair_num, fp);
    }
}

/* get the ls-fitted plane of each layer as in SCHNAaP: for rise calculation etc */
static void layer_lsfit(long num_bp, char *bseq, long **pair_num, long **seidx,
                        char **AtomName, double **xyz, FILE * fp)
{
    char str[BUF512];
    double odist, ppos[4], ppos0[4], z[4], z0[4], *adist, **bxyz;
    double zave[4], ddorg[4];
    long BASIZE, i, j, k, num, num_batom;
    long batom[NUM_BASE_ATOMS];

    BASIZE = pair_num[1][0] * NUM_BASE_ATOMS;
    adist = dvector(1, BASIZE);
    bxyz = dmatrix(1, BASIZE, 1, 3);

    print_sep(fp, '*', 76);
    fprintf(fp, "Least-squares fitted plane to all base atoms for each base layer."
            " Nx, Ny, Nz\nare the normal vector to the ls-fitted plane, and Ox, Oy,"
            " Oz are the point\n");
    fprintf(fp, "the plane passes through (geometric center of all the base atoms)."
            " The first\nnumber in [] is the standard deviation of the vertical"
            " distances of all base\n");
    fprintf(fp, "atoms from the plane, which gives a measure of planarity."
            " The second number\nis the projection of d(Ox, Oy, Oz) vector onto"
            " mean (Nx, Ny, Nz)\n");
    fprintf(fp, "     base                   Nx      Ny      Nz         Ox       Oy"
            "       Oz\n");
    for (i = 1; i <= num_bp; i++) {
        stepstr(i, 0, bseq, pair_num, str);
        num = 0;
        for (j = 1; j <= pair_num[i][0]; j++) {
            k = pair_num[i][j];
            cehs_base_atoms(AtomName, seidx[k][1], seidx[k][2], &num_batom, batom);
            for (k = 1; k <= num_batom; k++)
                cpxyz(xyz[batom[k]], bxyz[num + k]);
            num += num_batom;
        }
        if (i > 1) {  /* keep a copy of previous ones */
            cpxyz(ppos, ppos0);
            cpxyz(z, z0);
        }
        ls_plane(bxyz, num, z, ppos, &odist, adist);
        odist = std_dvector(adist, num);
        fprintf(fp, "%s [%5.3f : ", str, odist);
        if (i > 1) {
            sumxyz(z, z0, zave);
            ddxyz(ppos0, ppos, ddorg);
            vec_norm(zave);
            fprintf(fp, "%5.2f]", fabs(dot(zave, ddorg)));
        } else
            fprintf(fp, " ----]");
        fprintf(fp, "%8.4f%8.4f%8.4f%9.2f%9.2f%9.2f\n", z[1], z[2], z[3], ppos[1],
                ppos[2], ppos[3]);
    }

    free_dvector(adist, 1, BASIZE);
    free_dmatrix(bxyz, 1, BASIZE, 1, 3);
}

/* process a structure */
static void process_anyhelix(char *inpfile, long waters)
{
    char pdbfile[BUF512], outfile[BUF512];
    char *bseq, *ChainID, **AtomName, **ResName, **Miscs;
    double *chi, **orien, **org, **NC1xyz, **o3_p, **xyz;
    double *bp_orien, *bp_org;
    long i, j, eqnum, hetatm, num, num_bp, num_residue;
    long *RY, *ResSeq, **pair_num, **seidx, **idxCN;
    long *idx, **htm_water;
    FILE *fp;

    time_t run_time;

    /* get PDB file and pairing information */
    pair_num = read_pair(inpfile, pdbfile, outfile, &num_bp, &hetatm, &eqnum);

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

    /* get the numbering information of each residue */
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

    /* get base sequence */
    bseq = cvector(1, num_residue);
    RY = lvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);

    for (i = 1; i <= num_bp; i++)
        for (j = 1; j <= pair_num[i][0]; j++)
            if (pair_num[i][j] > num_residue || RY[pair_num[i][j]] < 0)
                fatal("residue number too big OR not a base!\n");

    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 3);
    NC1xyz = dmatrix(1, num_residue, 1, 7);  /* RN9/YN1 & C1' atomic coordinates */
    o3_p = dmatrix(1, num_residue, 1, 8);  /* O3'/P atomic coordinates */
    base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
              ResSeq, Miscs, xyz, orien, org, NC1xyz, o3_p);

    /* get chi torsion angle for each base residue */
    chi = dvector(1, num_residue);
    idxCN = lmatrix(1, num_residue, 1, 2);
    get_chi_angle(num_residue, RY, bseq, seidx, xyz, AtomName, ResName, ChainID,
                  ResSeq, Miscs, chi, idxCN);

    fp = open_file(outfile, "w");
    fprintf(fp, "    %s\n", Gvars.X3DNA_VER);
    print_sep(fp, '*', 76);
    fprintf(fp, "File name: %s\n", pdbfile);

    run_time = time(NULL);
    fprintf(fp, "Date and time: %s\n", ctime(&run_time));

    fprintf(fp, "Number of base-pairs: %ld\n", num_bp);
    fprintf(fp, "Number of atoms: %ld\n", num);

    print_sep(fp, '*', 76);
    stnd1_up(orien, org, chi, num_bp, pair_num, fp);

    print_sep(fp, '*', 76);
    print_pdb_title(pdbfile, "*", fp);
    print_sep(fp, '*', 76);

    bp_orien = dvector(1, num_bp * 9);
    bp_org = dvector(1, num_bp * 3);

    idx = lvector(1, num);
    atom_idx(num, AtomName, NULL, idx);

    htm_water = lmatrix(1, 4, 0, num);  /* HETATM and water index */
    init_htm_water(waters, num, num_residue, idx, htm_water);
    identify_htw(num_residue, seidx, RY, AtomName, ResName, ChainID, ResSeq,
                 Miscs, xyz, htm_water);

    cross_pair(num, num_residue, num_bp, pair_num, &Gvars.misc_pars, bseq, seidx,
               AtomName, ResName, ChainID, ResSeq, Miscs, xyz, chi, orien, org,
               NC1xyz, bp_orien, bp_org, RY, idx, htm_water, fp);
    print_ah_frame(num_bp, pair_num, bseq, bp_orien, bp_org);
    ah_step_hel(num_bp, pair_num, bseq, seidx, bp_orien, bp_org, AtomName, ResName,
                ChainID, ResSeq, Miscs, xyz, htm_water, fp);

    if (eqnum) {  /* all set has the same # of bases */
        local_c1par(num_bp, bseq, bp_orien, NC1xyz, pair_num, fp);
        get_ah_helix(num_bp, bseq, pair_num, num, idxCN, seidx, AtomName, xyz, NC1xyz, bp_org, fp);
    }

    layer_lsfit(num_bp, bseq, pair_num, seidx, AtomName, xyz, fp);

    close_file(fp);

    /* free allocated vectors & matrices */
    free_lmatrix(pair_num, 1, num_bp, 0, MBASES);
    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cvector(bseq, 1, num_residue);
    free_lvector(RY, 1, num_residue);
    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 3);
    free_dmatrix(NC1xyz, 1, num_residue, 1, 7);
    free_dmatrix(o3_p, 1, num_residue, 1, 8);
    free_dvector(chi, 1, num_residue);
    free_lmatrix(idxCN, 1, num_residue, 1, 2);
    free_dvector(bp_orien, 1, num_bp * 9);
    free_dvector(bp_org, 1, num_bp * 3);
    free_lvector(idx, 1, num);
    free_lmatrix(htm_water, 1, 4, 0, num);
}

int main(int argc, char *argv[])
{
    long i, j, waters = 0;
    time_t time0;

    time(&time0);

    /* clean up files with fixed name */
    remove_file(POC_FILE);
    remove_file(MUL_FILE);
    remove_file(ALLP_FILE);
    remove_file(STACK_FILE);
    remove_file(HSTACK_FILE);
    remove_file(REF_FILE);

    set_my_globals(argv[0]);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'W')
                waters = 1;
            else
                anyhelix_usage();
    }
    if (argc == i)
        process_anyhelix("stdin", waters);
    else
        for (j = i; j < argc; j++) {
            fprintf(stderr, "\n......Processing structure #%ld: <%s>......\n", j - i + 1, argv[j]);
            process_anyhelix(argv[j], waters);
        }

    clear_my_globals();

    print_used_time(time0);

    return 0;
}
