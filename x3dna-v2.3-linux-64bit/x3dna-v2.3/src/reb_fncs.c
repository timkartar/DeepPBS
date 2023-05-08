#include "x3dna.h"

/* link same strand O3'(i) to P (i+1) */
void link_o3_p(long num_residue, long **seidx, char **AtomName, double **xyz,
               char *ChainID, long **connect)
{
    double d_o3_p;
    long i, ib, ie, idx[7];
    long *O3, *P;

    O3 = lvector(1, num_residue);
    P = lvector(1, num_residue);

    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];
        O3[i] = find_1st_atom(" O3'", AtomName, ib, ie, "");
        P[i] = find_1st_atom(" P  ", AtomName, ib, ie, "");
    }

    /* connect same chain O3'(i) to P(i+1) if not yet */
    for (i = 1; i <= num_residue - 1; i++) {
        ib = O3[i];
        ie = P[i + 1];
        if (ib && ie && ChainID[ib] == ChainID[ie]) {
            d_o3_p = p1p2_dist(xyz[ib], xyz[ie]);
            if (d_o3_p < Gvars.misc_pars.o3p_dist) {
                if (connect[ib][7] + 1 > 6) {
                    fprintf(stderr, "atom <%ld: [%s]> has over 6 bonds\n", ib, AtomName[ib]);
                    continue;
                } else {
                    ++connect[ib][7];
                    connect[ib][connect[ib][7]] = ie;
                }
                lsort(connect[ib][7], connect[ib], idx);  /* re-order it */
                if (connect[ie][7] + 1 > 6) {
                    fprintf(stderr, "atom <%ld: [%s]> has over 6 bonds\n", ie, AtomName[ie]);
                    continue;
                } else {
                    ++connect[ie][7];
                    connect[ie][connect[ie][7]] = ib;
                }
                lsort(connect[ie][7], connect[ie], idx);  /* re-order it */
            } else
                fprintf(stderr, "O3' (#%ld) and P (#%ld) on chain %c have "
                        "distance %.1f over %.1f: no linkage assigned\n",
                        ib, ie, ChainID[ib], d_o3_p, Gvars.misc_pars.o3p_dist);
        }
    }

    free_lvector(O3, 1, num_residue);
    free_lvector(P, 1, num_residue);
}

/* get atom linkage information for rebuilt atomic model */
void atom_lkg(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
              double **xyz, long xml, char *outfile)
{
    long num_residue;
    long *RY, **connect, **seidx;

    if (xml) {  /* do NOT know how to add CONECT records in PDBML yet */
        write_pdbml(xml, num, AtomName, ResName, ChainID, ResSeq, xyz, outfile);
        return;
    }

    /* get the numbering information of each residue */
    seidx = residue_idx(num, ResSeq, NULL, ChainID, ResName, &num_residue);

    RY = lvector(1, num_residue);  /* default to 0: all nucleic acids */

    /* maximum 6 bonds, 7 holds #s, 8 holds the atom # in question */
    connect = lmatrix(1, num, 1, 7);
    get_bonds(num, AtomName, xyz, num_residue, RY, seidx, connect);

    link_o3_p(num_residue, seidx, AtomName, xyz, ChainID, connect);

    write_pdbcnt(num, AtomName, ResName, ChainID, ResSeq, xyz, connect, outfile);

    free_lvector(RY, 1, num_residue);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_lmatrix(connect, 1, num, 1, 7);
}

/* build a single strand DNA structure in PDB format */
void atomic_pdb1(long num_bp, long num_atoms, long num_max_per_residue, long is_helical,
                 long xdir, char **bp_seq, double **step_par, char *BDIR, long xml, char *outfile)
{
    char spdb[BUF512];
    char *sChainID, *tChainID;
    char **sAtomName, **sResName, **tAtomName, **tResName;
    double mpos[4], pos[4];
    double pos_next[4] = { EMPTY_NUMBER, 0.0, 0.0, 0.0 };
    double **mst, **orien, **orien_next, **sxyz, **txyz;
    long i, ik, j, k, snum, tnum = 0;
    long *sResSeq, *tResSeq;
    FILE *fp;

    tAtomName = cmatrix(1, num_atoms, 0, 4);
    tResName = cmatrix(1, num_atoms, 0, 3);
    tChainID = cvector(1, num_atoms);
    tResSeq = lvector(1, num_atoms);
    txyz = dmatrix(1, num_atoms, 1, 3);

    sAtomName = cmatrix(1, num_max_per_residue, 0, 4);
    sResName = cmatrix(1, num_max_per_residue, 0, 3);
    sChainID = cvector(1, num_max_per_residue);
    sResSeq = lvector(1, num_max_per_residue);
    sxyz = dmatrix(1, num_max_per_residue, 1, 3);

    mst = dmatrix(1, 3, 1, 3);
    orien = dmatrix(1, 3, 1, 3);
    orien_next = dmatrix(1, 3, 1, 3);
    identity_matrix(orien_next, 3);

    fp = open_file(REF_FILE, "w");
    fprintf(fp, "%5ld bases\n", num_bp);

    for (i = 1; i <= num_bp; i++) {
        set_std_base_pdb(BDIR, FALSE, bp_seq[i][1], spdb);
        snum = read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL, 1, "*");
        if (xdir)
            for (j = 1; j <= snum; j++) {
                sxyz[j][1] = -sxyz[j][1];  /* reverse x */
                sxyz[j][3] = -sxyz[j][3];  /* reverse z */
            }
        if (is_helical)
            xhelfunc(step_par[i], orien, mst, pos, mpos);
        else
            xbpfunc(step_par[i], orien, mst, pos, mpos);

        multi_vec_Tmatrix(pos, 3, orien_next, 3, 3, mpos);
        sumxyz(mpos, pos_next, pos_next);

        multi_matrix(orien_next, 3, 3, orien, 3, 3, mst);
        copy_dmatrix(mst, 3, 3, orien_next);

        for (j = 1; j <= snum; j++) {
            ik = tnum + j;
            strcpy(tAtomName[ik], sAtomName[j]);
            sprintf(tResName[ik], "  %c", bp_seq[i][1]);
            tChainID[ik] = Gvars.REBUILD_CHAIN_IDS[0];
            tResSeq[ik] = i;
            for (k = 1; k <= 3; k++)
                txyz[ik][k] = dot(sxyz[j], orien_next[k]) + pos_next[k];
        }

        tnum += snum;

        /* write out reference frames */
        fprintf(fp, "... %5ld %c ...\n", i, bp_seq[i][1]);
        print_ref_frames(fp, pos_next, orien_next);
    }

    close_file(fp);

    atom_lkg(tnum, tAtomName, tResName, tChainID, tResSeq, txyz, xml, outfile);

    free_pdb(num_atoms, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
    free_pdb(num_max_per_residue, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(orien, 1, 3, 1, 3);
    free_dmatrix(orien_next, 1, 3, 1, 3);
}

/* build a duplex DNA structure in PDB format */
void atomic_pdb2(long parallel, long num_bp, long num_atoms, long num_max_per_residue,
                 long is_helical, long xdir, char **bp_seq, double **bp_par,
                 double **step_par, char *BDIR, long xml, char *outfile)
{
    char fname1[BUF512], fname2[BUF512];
    char *tChainID;
    char **sAtomName1, **sAtomName2, **tAtomName, **tAtomName2, **tResName;
    double mpos[4], pos[4];
    double pos_next[4] = { EMPTY_NUMBER, 0.0, 0.0, 0.0 };
    double **mst, **orien, **orien_next, **sxyz1, **sxyz2, **txyz, **txyz2;
    long i, ik, j, k, num1, num2, tnum = 0, tnum2 = 0;
    long *tResSeq;
    long **s2idx;  /* strand II index */
    FILE *fp;

    tAtomName = cmatrix(1, num_atoms, 0, 4);
    tResName = cmatrix(1, num_atoms, 0, 3);
    tChainID = cvector(1, num_atoms);
    tResSeq = lvector(1, num_atoms);
    txyz = dmatrix(1, num_atoms, 1, 3);
    tAtomName2 = cmatrix(1, num_atoms, 0, 4);
    txyz2 = dmatrix(1, num_atoms, 1, 3);
    sAtomName1 = cmatrix(1, num_max_per_residue, 0, 4);
    sxyz1 = dmatrix(1, num_max_per_residue, 1, 3);
    sAtomName2 = cmatrix(1, num_max_per_residue, 0, 4);
    sxyz2 = dmatrix(1, num_max_per_residue, 1, 3);

    mst = dmatrix(1, 3, 1, 3);
    orien = dmatrix(1, 3, 1, 3);
    orien_next = dmatrix(1, 3, 1, 3);
    identity_matrix(orien_next, 3);

    s2idx = lmatrix(1, num_bp, 1, 2);

    fp = open_file(REF_FILE, "w");
    fprintf(fp, "%5ld base-pairs\n", num_bp);

    for (i = 1; i <= num_bp; i++) {
        set_std_base_pdb(BDIR, FALSE, bp_seq[i][1], fname1);
        set_std_base_pdb(BDIR, FALSE, bp_seq[i][2], fname2);
        set_bp_pdb(num_max_per_residue, fname1, fname2, bp_par[i], xdir, &num1, &num2,
                   sAtomName1, sAtomName2, sxyz1, sxyz2, bp_seq[i][0]);

        if (is_helical)
            xhelfunc(step_par[i], orien, mst, pos, mpos);
        else
            xbpfunc(step_par[i], orien, mst, pos, mpos);

        multi_vec_Tmatrix(pos, 3, orien_next, 3, 3, mpos);
        sumxyz(mpos, pos_next, pos_next);

        multi_matrix(orien_next, 3, 3, orien, 3, 3, mst);
        copy_dmatrix(mst, 3, 3, orien_next);

        for (j = 1; j <= num1; j++) {
            ik = tnum + j;
            strcpy(tAtomName[ik], sAtomName1[j]);
            sprintf(tResName[ik], "  %c", bp_seq[i][1]);
            tChainID[ik] = Gvars.REBUILD_CHAIN_IDS[0];
            tResSeq[ik] = i;
            for (k = 1; k <= 3; k++)
                txyz[ik][k] = dot(sxyz1[j], orien_next[k]) + pos_next[k];
        }

        for (j = 1; j <= num2; j++) {
            ik = tnum2 + j;
            strcpy(tAtomName2[ik], sAtomName2[j]);
            for (k = 1; k <= 3; k++)
                txyz2[ik][k] = dot(sxyz2[j], orien_next[k]) + pos_next[k];
        }

        s2idx[i][1] = tnum2;
        s2idx[i][2] = num2;

        tnum += num1;
        tnum2 += num2;

        /* write out reference frames */
        fprintf(fp, "... %5ld %c%c%c ...\n", i, bp_seq[i][1], bp_seq[i][0], bp_seq[i][2]);
        print_ref_frames(fp, pos_next, orien_next);
    }

    close_file(fp);

    if (parallel)
        combine_pstnd2(num_bp, bp_seq, s2idx, &tnum, tAtomName, tResName,
                       tChainID, tResSeq, txyz, tAtomName2, txyz2);
    else  /* anti-parallel: "reverse" strand II and combined with I */
        reverse_stnd2(num_bp, bp_seq, s2idx, &tnum, tAtomName, tResName,
                      tChainID, tResSeq, txyz, tAtomName2, txyz2, 0);

    atom_lkg(tnum, tAtomName, tResName, tChainID, tResSeq, txyz, xml, outfile);

    free_pdb(num_atoms, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
    free_cmatrix(tAtomName2, 1, num_atoms, 0, 4);
    free_dmatrix(txyz2, 1, num_atoms, 1, 3);
    free_cmatrix(sAtomName1, 1, num_max_per_residue, 0, 4);
    free_dmatrix(sxyz1, 1, num_max_per_residue, 1, 3);
    free_cmatrix(sAtomName2, 1, num_max_per_residue, 0, 4);
    free_dmatrix(sxyz2, 1, num_max_per_residue, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(orien, 1, 3, 1, 3);
    free_dmatrix(orien_next, 1, 3, 1, 3);
    free_lmatrix(s2idx, 1, num_bp, 1, 2);
}

/* get base & C1' atoms index in a residue based on atom name */
void base_c1_atoms(char **AtomName, long ib, long ie, long *num_batom, long *batom)
{
    long i;

    *num_batom = 0;

    for (i = ib; i <= ie; i++)
        if (is_baseatom(AtomName[i]) || (!strcmp(AtomName[i], " C1'"))) {
            ++*num_batom;
            if (*num_batom > NUM_BASE_ATOMS)
                fatal("too many base atoms in a residue\n");
            batom[*num_batom] = i;
        }
}

void extract_base_atoms(long num, char **AtomName, double **xyz, long *bnum,
                        char **bAtomName, double **bxyz)
{
    long i, *batom;

    batom = lvector(1, NUM_BASE_ATOMS);

    base_c1_atoms(AtomName, 1, num, bnum, batom);
    for (i = 1; i <= *bnum; i++) {
        strcpy(bAtomName[i], AtomName[batom[i]]);
        cpxyz(xyz[batom[i]], bxyz[i]);
    }

    free_lvector(batom, 1, NUM_BASE_ATOMS);
}

/* build a DNA structure with only base and P atoms. if the residue contains
   sugar-phosphate backbone, they will be deleted. thus the generated
   structure will have less than [num_atoms + 2 * (num_bp - 1)] atoms */
void atomic_base_p(long parallel, long num_bp, long num_atoms, long num_max_per_residue,
                   long is_helical, long xdir, char **bp_seq, double **bp_par,
                   double **step_par, char *BDIR, long *pidx, long xml, char *outfile)
{
    char *Patom = " P  ";
    char fname1[BUF512], fname2[BUF512], str[BUF512];
    char *tChainID;
    char **bAtomName1, **bAtomName2, **sAtomName1, **sAtomName2;
    char **tAtomName, **tAtomName2, **tResName;
    double mpos[4], mst_pos[4], pos[4];
    double pos_next[4] = { EMPTY_NUMBER, 0.0, 0.0, 0.0 };
    double **mst, **mst_next, **orien, **orien_next;
    double **bxyz1, **bxyz2, **sxyz1, **sxyz2;
    double **ppos, **txyz, **txyz2;
    long bnum1, bnum2, snum1, snum2;
    long i, ik, j, k, tnum = 0, tnum2 = 0;
    long num_pidx = 2048;  /* assumed distinct P xyz index */
    long *tResSeq;
    long **s2idx;  /* strand II index */
    FILE *fp;

    ppos = dmatrix(1, num_pidx, 1, 6);
    init_dmatrix(ppos, 1, num_pidx, 1, 6, EMPTY_NUMBER);

    sprintf(fname1, "Pxyz%s.dat", (is_helical) ? "H" : "");
    get_BDIR(str, fname1);
    strcat(str, fname1);

    fp = open_file(str, "r");
    fprintf(stderr, " ...... reading file: %s ...... \n", str);

    for (i = 1; i <= 6; i++)  /* skip 6 title lines */
        if (fgets(str, sizeof str, fp) == NULL)
            fatal("error reading parameter file: six title lines\n");

    while (fgets(str, sizeof str, fp) != NULL) {
        if (sscanf(str, "%ld", &i) != 1)
            fatal("error reading P xyz index\n");
        if (!lval_in_range(i, 1, num_pidx))
            fatal("P xyz index out of range\n");
        ik = sscanf(str, "%*d %lf %lf %lf %lf %lf %lf", &ppos[i][1],
                    &ppos[i][2], &ppos[i][3], &ppos[i][4], &ppos[i][5], &ppos[i][6]);
        if (ik == 3) {
            if (parallel) {  /* e.g., Hoogsteen: not symmetrical */
                fprintf(stderr, "Only 3 coordinates provided. Copied used for the" " other 3\n");
                ppos[i][4] = ppos[i][1];
                ppos[i][5] = ppos[i][2];
                ppos[i][6] = ppos[i][3];
            } else {  /* anti-parallel: reverse y and z */
                ppos[i][4] = ppos[i][1];
                ppos[i][5] = -ppos[i][2];
                ppos[i][6] = -ppos[i][3];
            }
        } else if (ik == 6);  /* all six, do nothing */
        else
            fatal("only 3 or 6 P-xyz-mst coordinates allowed per line\n");
    }

    close_file(fp);

    num_atoms += 2 * (num_bp - 1);  /* including P atoms */

    tAtomName = cmatrix(1, num_atoms, 0, 4);
    tResName = cmatrix(1, num_atoms, 0, 3);
    tChainID = cvector(1, num_atoms);
    tResSeq = lvector(1, num_atoms);
    txyz = dmatrix(1, num_atoms, 1, 3);
    tAtomName2 = cmatrix(1, num_atoms, 0, 4);
    txyz2 = dmatrix(1, num_atoms, 1, 3);
    sAtomName1 = cmatrix(1, num_max_per_residue, 0, 4);
    sxyz1 = dmatrix(1, num_max_per_residue, 1, 3);
    sAtomName2 = cmatrix(1, num_max_per_residue, 0, 4);
    sxyz2 = dmatrix(1, num_max_per_residue, 1, 3);
    bAtomName1 = cmatrix(1, NUM_BASE_ATOMS, 0, 4);
    bxyz1 = dmatrix(1, NUM_BASE_ATOMS, 1, 3);
    bAtomName2 = cmatrix(1, NUM_BASE_ATOMS, 0, 4);
    bxyz2 = dmatrix(1, NUM_BASE_ATOMS, 1, 3);

    mst = dmatrix(1, 3, 1, 3);
    mst_next = dmatrix(1, 3, 1, 3);
    identity_matrix(mst_next, 3);

    orien = dmatrix(1, 3, 1, 3);
    orien_next = dmatrix(1, 3, 1, 3);
    identity_matrix(orien_next, 3);

    s2idx = lmatrix(1, num_bp, 1, 2);

    fp = open_file(REF_FILE, "w");
    fprintf(fp, "%5ld base-pairs\n", num_bp);

    for (i = 1; i <= num_bp; i++) {
        set_std_base_pdb(BDIR, FALSE, bp_seq[i][1], fname1);
        set_std_base_pdb(BDIR, FALSE, bp_seq[i][2], fname2);
        set_bp_pdb(num_max_per_residue, fname1, fname2, bp_par[i], xdir, &snum1,
                   &snum2, sAtomName1, sAtomName2, sxyz1, sxyz2, bp_seq[i][0]);

        extract_base_atoms(snum1, sAtomName1, sxyz1, &bnum1, bAtomName1, bxyz1);
        extract_base_atoms(snum2, sAtomName2, sxyz2, &bnum2, bAtomName2, bxyz2);

        if (i == 1 && snum1 > bnum1 && snum2 > bnum2)
            fprintf(stderr, "ignoring sugar-phosphate atoms in your bases\n");

        if (is_helical)
            xhelfunc(step_par[i], orien, mst, pos, mpos);
        else
            xbpfunc(step_par[i], orien, mst, pos, mpos);

        if (i > 1) {  /* attach P atoms */
            if (ppos[pidx[i]][1] < EMPTY_CRITERION) {
                fprintf(stderr, "P xyz index %ld does not exist."
                        " use default for B-DNA instead\n", pidx[i]);
                pidx[i] = 2;
            }
            multi_vec_Tmatrix(mpos, 3, orien_next, 3, 3, mst_pos);
            sumxyz(pos_next, mst_pos, mst_pos);

            multi_matrix(orien_next, 3, 3, mst, 3, 3, mst_next);

            tnum++;
            strcpy(tAtomName[tnum], Patom);
            sprintf(tResName[tnum], "  %c", bp_seq[i][1]);
            tChainID[tnum] = 'A';
            tResSeq[tnum] = i;
            for (j = 1; j <= 3; j++)
                txyz[tnum][j] = dot(ppos[pidx[i]], mst_next[j]) + mst_pos[j];

            ik = tnum2 + bnum2 + 1;
            strcpy(tAtomName2[ik], Patom);
            for (j = 1; j <= 3; j++)
                txyz2[ik][j] = dot(ppos[pidx[i]] + 3, mst_next[j]) + mst_pos[j];
        }
        multi_vec_Tmatrix(pos, 3, orien_next, 3, 3, mpos);
        sumxyz(mpos, pos_next, pos_next);

        multi_matrix(orien_next, 3, 3, orien, 3, 3, mst);
        copy_dmatrix(mst, 3, 3, orien_next);

        for (j = 1; j <= bnum1; j++) {
            ik = tnum + j;
            strcpy(tAtomName[ik], bAtomName1[j]);
            sprintf(tResName[ik], "  %c", bp_seq[i][1]);
            tChainID[ik] = 'A';
            tResSeq[ik] = i;
            for (k = 1; k <= 3; k++)
                txyz[ik][k] = dot(bxyz1[j], orien_next[k]) + pos_next[k];
        }

        for (j = 1; j <= bnum2; j++) {
            ik = tnum2 + j;
            strcpy(tAtomName2[ik], bAtomName2[j]);
            for (k = 1; k <= 3; k++)
                txyz2[ik][k] = dot(bxyz2[j], orien_next[k]) + pos_next[k];
        }

        s2idx[i][1] = tnum2;
        s2idx[i][2] = (i > 1) ? bnum2 + 1 : bnum2;

        tnum += bnum1;
        tnum2 += (i > 1) ? bnum2 + 1 : bnum2;

        /* write out reference frames */
        fprintf(fp, "... %5ld %c%c%c ...\n", i, bp_seq[i][1], bp_seq[i][0], bp_seq[i][2]);
        print_ref_frames(fp, pos_next, orien_next);
    }

    close_file(fp);

    if (parallel)  /* very UNLIKELY to be used! */
        combine_pstnd2(num_bp, bp_seq, s2idx, &tnum, tAtomName, tResName,
                       tChainID, tResSeq, txyz, tAtomName2, txyz2);
    else  /* anti-parallel: "reverse" strand II and combined with I */
        reverse_stnd2(num_bp, bp_seq, s2idx, &tnum, tAtomName, tResName,
                      tChainID, tResSeq, txyz, tAtomName2, txyz2, 1);
    atom_lkg(tnum, tAtomName, tResName, tChainID, tResSeq, txyz, xml, outfile);

    free_dmatrix(ppos, 1, num_pidx, 1, 6);
    free_pdb(num_atoms, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
    free_cmatrix(tAtomName2, 1, num_atoms, 0, 4);
    free_dmatrix(txyz2, 1, num_atoms, 1, 3);
    free_cmatrix(sAtomName1, 1, num_max_per_residue, 0, 4);
    free_dmatrix(sxyz1, 1, num_max_per_residue, 1, 3);
    free_cmatrix(sAtomName2, 1, num_max_per_residue, 0, 4);
    free_dmatrix(sxyz2, 1, num_max_per_residue, 1, 3);
    free_cmatrix(bAtomName1, 1, NUM_BASE_ATOMS, 0, 4);
    free_dmatrix(bxyz1, 1, NUM_BASE_ATOMS, 1, 3);
    free_cmatrix(bAtomName2, 1, NUM_BASE_ATOMS, 0, 4);
    free_dmatrix(bxyz2, 1, NUM_BASE_ATOMS, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(mst_next, 1, 3, 1, 3);
    free_dmatrix(orien, 1, 3, 1, 3);
    free_dmatrix(orien_next, 1, 3, 1, 3);
    free_lmatrix(s2idx, 1, num_bp, 1, 2);
}

/* get total number of atoms in the rebuilt structure for allocation */
void num_PDB_atoms(long num_bp, long is_single, char **bp_seq, char *BDIR,
                   long *num_atoms, long *num_max_per_residue)
{
    char spdb[BUF512];
    long kmax = -10000, kmin = 10000;
    long ds, i, j, k, tnum = 0, pmax = -10000;
    long *idx;

    ds = (is_single) ? 1 : 2;

    for (i = 1; i <= num_bp; i++)
        for (j = 1; j <= ds; j++) {
            if (kmax < bp_seq[i][j])
                kmax = bp_seq[i][j];
            if (kmin > bp_seq[i][j])
                kmin = bp_seq[i][j];
        }
    if (kmax < 0 || kmin < 1)
        fatal("there are wrong base types\n");

    idx = lvector(1, kmax);
    for (i = 1; i <= kmax; i++)
        idx[i] = 0;

    for (i = 1; i <= num_bp; i++)
        for (j = 1; j <= ds; j++) {
            ++idx[(long) bp_seq[i][j]];
        }

    for (i = 1; i <= kmax; i++)
        if (idx[i]) {
            set_std_base_pdb(BDIR, FALSE, (char) i, spdb);
            k = number_of_atoms(spdb, 1, "*");
            if (k > pmax)
                pmax = k;
            tnum += idx[i] * k;
        }
    free_lvector(idx, 1, kmax);

    *num_atoms = tnum;
    *num_max_per_residue = pmax;
}

/* set a base-pair (PDB format) w.r.t. its reference frame */
void set_bp_pdb(long num_max_per_residue, char *fname1, char *fname2, double *param,
                long xdir, long *num1, long *num2, char **AtomName1, char **AtomName2,
                double **xyz1, double **xyz2, char ap)
{
    char *ChainID, **ResName;
    double mpos[4], pos[4], **mst, **temp, **orien;
    long i, j, *ResSeq;

    ResName = cmatrix(1, num_max_per_residue, 0, 3);
    ChainID = cvector(1, num_max_per_residue);
    ResSeq = lvector(1, num_max_per_residue);

    *num1 = read_pdb(fname1, NULL, AtomName1, ResName, ChainID, ResSeq, xyz1, NULL, 1, "*");
    *num2 = read_pdb(fname2, NULL, AtomName2, ResName, ChainID, ResSeq, xyz2, NULL, 1, "*");

    if (xdir) {  /* reverse x- and z-axes */
        for (i = 1; i <= *num1; i++) {
            xyz1[i][1] = -xyz1[i][1];
            xyz1[i][3] = -xyz1[i][3];
        }
        for (i = 1; i <= *num2; i++) {
            xyz2[i][1] = -xyz2[i][1];
            xyz2[i][3] = -xyz2[i][3];
        }
    }
    if (ap == '-')  /* anti-parallel: reverse y- and z-axes for base II */
        for (i = 1; i <= *num2; i++) {
            xyz2[i][2] = -xyz2[i][2];
            xyz2[i][3] = -xyz2[i][3];
        }

    mst = dmatrix(1, 3, 1, 3);
    orien = dmatrix(1, 3, 1, 3);
    xbpfunc(param, orien, mst, pos, mpos);

    temp = dmatrix(1, num_max_per_residue, 1, 3);

    for (i = 1; i <= *num1; i++)
        for (j = 1; j <= 3; j++)
            temp[i][j] = dot(xyz1[i], orien[j]) + pos[j];
    copy_dmatrix(temp, *num1, 3, xyz1);
    change_xyz(0, mpos, mst, *num1, xyz1);
    change_xyz(0, mpos, mst, *num2, xyz2);

    free_pdb(num_max_per_residue, NULL, NULL, ResName, ChainID, ResSeq, NULL, NULL);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(orien, 1, 3, 1, 3);
    free_dmatrix(temp, 1, num_max_per_residue, 1, 3);
}

/* get R or Y base file name */
void base_fname(char bname, char *BDIR, char *fname)
{
    static char *Rbases = "AGI";

    if (strchr(Rbases, toupper((int) bname)) != NULL)
        sprintf(fname, "%sBlock_R.alc", BDIR);
    else
        sprintf(fname, "%sBlock_Y.alc", BDIR);
}

/* get schematic representation of DNA in ALCHEMY format */
void block_alc1(long num_bp, long is_single, long is_helical, long xdir, char **bp_seq,
                double **step_par, char *BDIR, char *outfile)
{
    char bpi[3], fname[BUF512];
    char **AtomName, **tAtomName;
    double mpos[4], pos[4];
    double pos_next[4] = { EMPTY_NUMBER, 0.0, 0.0, 0.0 };
    double **mst, **org_xyz, **orien, **orien_next, **txyz, **xyz;
    long i, ik, ia = 0, ib = 0, j, k;
    long bidx, nbond, num, tnbond, tnum;
    long *ibase, *tibase;
    long **linkage, **tlinkage;
    FILE *fp;

    strcpy(fname, BDIR);
    strcat(fname, "Block_R.alc");
    get_alc_nums(fname, &num, &nbond);

    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 1, 2);

    tnum = (num + 1) * num_bp;
    tnbond = (nbond + 1) * num_bp - 1;
    tAtomName = cmatrix(1, tnum, 0, 2);
    txyz = dmatrix(1, tnum, 1, 3);
    tibase = lvector(1, tnum);
    tlinkage = lmatrix(1, tnbond, 1, 2);

    mst = dmatrix(1, 3, 1, 3);
    orien = dmatrix(1, 3, 1, 3);
    orien_next = dmatrix(1, 3, 1, 3);
    identity_matrix(orien_next, 3);

    org_xyz = dmatrix(1, num_bp, 1, 3);

    fp = open_file(REF_FILE, "w");
    if (is_single)
        fprintf(fp, "%5ld bases\n", num_bp);
    else
        fprintf(fp, "%5ld base-pairs\n", num_bp);

    for (i = 1; i <= num_bp; i++) {
        if (is_single)
            base_idx(i, &bp_seq[i][1], &bidx, 1);
        else {
            sprintf(bpi, "%c%c", toupper((int) bp_seq[i][1]), toupper((int) bp_seq[i][2]));
            bidx = basepair_idx(bpi);
        }

        if (is_single)
            base_fname(bp_seq[i][1], BDIR, fname);
        else
            sprintf(fname, "%sBlock_BP.alc", BDIR);

        read_alc(fname, &num, &nbond, AtomName, xyz, ibase, linkage);

        if (xdir)
            for (j = 1; j <= num; j++) {
                xyz[j][1] = -xyz[j][1];  /* reverse x */
                xyz[j][3] = -xyz[j][3];  /* reverse z */
            }
        if (is_helical)
            xhelfunc(step_par[i], orien, mst, pos, mpos);
        else
            xbpfunc(step_par[i], orien, mst, pos, mpos);

        multi_vec_Tmatrix(pos, 3, orien_next, 3, 3, mpos);
        sumxyz(mpos, pos_next, pos_next);

        multi_matrix(orien_next, 3, 3, orien, 3, 3, mst);
        copy_dmatrix(mst, 3, 3, orien_next);

        for (j = 1; j <= num; j++) {
            ik = ia + j;
            strcpy(tAtomName[ik], AtomName[j]);
            tibase[ik] = bidx;
            for (k = 1; k <= 3; k++)
                txyz[ik][k] = dot(xyz[j], orien_next[k]) + pos_next[k];
        }

        for (j = 1; j <= nbond; j++) {
            ik = ib + j;
            for (k = 1; k <= 2; k++)
                tlinkage[ik][k] = ia + linkage[j][k];
        }

        ia += num;
        ib += nbond;

        cpxyz(pos_next, org_xyz[i]);

        /* write out reference frames */
        if (is_single)
            fprintf(fp, "... %5ld %c ...\n", i, bp_seq[i][1]);
        else
            fprintf(fp, "... %5ld %c%c%c ...\n", i, bp_seq[i][1], bp_seq[i][0], bp_seq[i][2]);
        print_ref_frames(fp, pos_next, orien_next);
    }

    close_file(fp);

    cnct_org(num_bp, ia, ib, tAtomName, txyz, tibase, tlinkage, org_xyz);

    write_alc(tnum, tnbond, tAtomName, txyz, tibase, tlinkage, outfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_alc(tnum, tnbond, tAtomName, txyz, tibase, 1, tlinkage);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(orien, 1, 3, 1, 3);
    free_dmatrix(orien_next, 1, 3, 1, 3);
    free_dmatrix(org_xyz, 1, num_bp, 1, 3);
}

void block_alc2(long num_bp, long is_helical, long xdir, char **bp_seq, double **bp_par,
                double **step_par, char *BDIR, char *outfile)
/* get schematic representation of a duplex DNA structure in ALCHEMY format.
   two blocks per base-pair, R & Y could have the same or different sizes */
{
    char fname1[BUF512], fname2[BUF512];
    char **AtomName, **tAtomName;
    double mpos[4], pos[4];
    double pos_next[4] = { EMPTY_NUMBER, 0.0, 0.0, 0.0 };
    double **mst, **orien, **orien_next;
    double **bp_xyz, **org_xyz, **txyz, **xyz;
    long i, ik, ia = 0, ib = 0, j, k;
    long bidx1, bidx2, bp_num, bp_nbond, nbond, num, tnbond, tnum;
    long *ibase, *tibase;
    long **linkage, **tlinkage;
    FILE *fp;

    strcpy(fname1, BDIR);
    strcat(fname1, "Block_R.alc");
    get_alc_nums(fname1, &num, &nbond);
    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 1, 2);
    read_alc(fname1, &num, &nbond, AtomName, xyz, ibase, linkage);

    bp_num = 2 * num;
    bp_nbond = 2 * nbond;
    tnum = (bp_num + 1) * num_bp;
    tnbond = (bp_nbond + 1) * num_bp - 1;

    bp_xyz = dmatrix(1, bp_num, 1, 3);

    tAtomName = cmatrix(1, tnum, 0, 2);
    txyz = dmatrix(1, tnum, 1, 3);
    tibase = lvector(1, tnum);
    tlinkage = lmatrix(1, tnbond, 1, 2);

    mst = dmatrix(1, 3, 1, 3);
    orien = dmatrix(1, 3, 1, 3);
    orien_next = dmatrix(1, 3, 1, 3);
    identity_matrix(orien_next, 3);

    org_xyz = dmatrix(1, num_bp, 1, 3);

    fp = open_file(REF_FILE, "w");
    fprintf(fp, "%5ld base-pairs\n", num_bp);

    for (i = 1; i <= num_bp; i++) {
        base_idx(i, &bp_seq[i][1], &bidx1, 1);
        base_idx(i, &bp_seq[i][2], &bidx2, 1);
        base_fname(bp_seq[i][1], BDIR, fname1);
        base_fname(bp_seq[i][2], BDIR, fname2);

        set_bp_alc(fname1, fname2, bp_par[i], xdir, num, nbond, bp_xyz, bp_seq[i][0]);

        (is_helical) ? xhelfunc(step_par[i], orien, mst, pos, mpos) :
            xbpfunc(step_par[i], orien, mst, pos, mpos);

        multi_vec_Tmatrix(pos, 3, orien_next, 3, 3, mpos);
        sumxyz(mpos, pos_next, pos_next);

        multi_matrix(orien_next, 3, 3, orien, 3, 3, mst);
        copy_dmatrix(mst, 3, 3, orien_next);

        for (j = 1; j <= num; j++) {  /* I */
            ik = ia + j;
            strcpy(tAtomName[ik], AtomName[j]);
            tibase[ik] = bidx1;
        }
        for (j = 1; j <= num; j++) {  /* II */
            ik = ia + num + j;
            strcpy(tAtomName[ik], AtomName[j]);
            tibase[ik] = bidx2;
        }

        for (j = 1; j <= bp_num; j++) {
            ik = ia + j;
            for (k = 1; k <= 3; k++)
                txyz[ik][k] = dot(bp_xyz[j], orien_next[k]) + pos_next[k];
        }

        for (j = 1; j <= nbond; j++) {  /* I */
            ik = ib + j;
            for (k = 1; k <= 2; k++)
                tlinkage[ik][k] = ia + linkage[j][k];
        }
        for (j = 1; j <= nbond; j++) {  /* II */
            ik = ib + nbond + j;
            for (k = 1; k <= 2; k++)
                tlinkage[ik][k] = ia + num + linkage[j][k];
        }

        ia += bp_num;
        ib += bp_nbond;

        cpxyz(pos_next, org_xyz[i]);

        /* write out reference frames */
        fprintf(fp, "... %5ld %c%c%c ...\n", i, bp_seq[i][1], bp_seq[i][0], bp_seq[i][2]);
        print_ref_frames(fp, pos_next, orien_next);
    }

    close_file(fp);

    cnct_org(num_bp, ia, ib, tAtomName, txyz, tibase, tlinkage, org_xyz);

    write_alc(tnum, tnbond, tAtomName, txyz, tibase, tlinkage, outfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_dmatrix(bp_xyz, 1, bp_num, 1, 3);
    free_alc(tnum, tnbond, tAtomName, txyz, tibase, 1, tlinkage);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(orien, 1, 3, 1, 3);
    free_dmatrix(orien_next, 1, 3, 1, 3);
    free_dmatrix(org_xyz, 1, num_bp, 1, 3);
}

/* set a base-pair (ALCHEMY format) w.r.t. its reference frame */
void set_bp_alc(char *fname1, char *fname2, double *param, long xdir, long num,
                long nbond, double **bp_xyz, char ap)
{
    char **AtomName;
    double mpos[4], pos[4];
    double **mst, **orien, **temp, **xyz2;
    long i, j;
    long *ibase;
    long **linkage;

    AtomName = cmatrix(1, num, 0, 2);
    temp = dmatrix(1, num, 1, 3);
    xyz2 = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 1, 2);

    read_alc(fname1, &num, &nbond, AtomName, bp_xyz, ibase, linkage);
    read_alc(fname2, &num, &nbond, AtomName, xyz2, ibase, linkage);

    if (xdir)  /* reverse x- and z-axes */
        for (i = 1; i <= num; i++) {
            bp_xyz[i][1] = -bp_xyz[i][1];
            bp_xyz[i][3] = -bp_xyz[i][3];
            xyz2[i][1] = -xyz2[i][1];
            xyz2[i][3] = -xyz2[i][3];
        }
    if (ap == '-')  /* anti-parallel: reverse y- and z-axes for base II */
        for (i = 1; i <= num; i++) {
            xyz2[i][2] = -xyz2[i][2];
            xyz2[i][3] = -xyz2[i][3];
        }

    mst = dmatrix(1, 3, 1, 3);
    orien = dmatrix(1, 3, 1, 3);
    xbpfunc(param, orien, mst, pos, mpos);

    for (i = 1; i <= num; i++)
        for (j = 1; j <= 3; j++)
            temp[i][j] = dot(bp_xyz[i], orien[j]) + pos[j];
    copy_dmatrix(temp, num, 3, bp_xyz);
    change_xyz(0, mpos, mst, num, bp_xyz);
    change_xyz(0, mpos, mst, num, xyz2);

    for (i = 1; i <= num; i++)
        cpxyz(xyz2[i], bp_xyz[i + num]);

    free_alc(num, nbond, AtomName, xyz2, ibase, 1, linkage);
    free_dmatrix(temp, 1, num, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(orien, 1, 3, 1, 3);
}

/* get the orientation and position of one base (base-pair) relative
   to the other given 6 base-pair or step parameters. the middle frame
   is also defined
   % ---------------------------------------------------------------
   %   param: 6 base-pair or step parameters in the order of:
   %            x       y        z       x        y         z
   %          Shear  Stretch  Stagger  Buckle Propeller  Opening
   %          Shift   Slide     Rise    Tilt    Roll      Twist
   % ---------------------------------------------------------------
   %   orien: base 2 (base-pair 2) orientation w.r.t. 1
   %   mst: middle frame orientation w.r.t. 1
   %   pos: base 2 (base-pair 2) origin w.r.t. 1
   %   mpos: middle frame origin w.r.t. 1
   % --------------------------------------------------------------- */
void xbpfunc(double *param, double **orien, double **mst, double *pos, double *mpos)
{
    double phi, yx_ang;
    double hinge[4];
    double y[4] = { EMPTY_NUMBER, 0.0, 1.0, 0.0 };
    double z[4] = { EMPTY_NUMBER, 0.0, 0.0, 1.0 };
    double **rmtx, **temp;
    long i;

    /* get the yx_ang angle (POSITIVE) and phi angle (-180, +180) */
    hinge[1] = param[4];
    hinge[2] = param[5];
    hinge[3] = 0.0;
    yx_ang = veclen(hinge);
    phi = vec_ang(hinge, y, z);

    rmtx = dmatrix(1, 3, 1, 3);
    temp = dmatrix(1, 3, 1, 3);

    /* base 2 (base-pair 2) orientation w.r.t. 1 */
    rotz(0.5 * param[6] - phi, orien);
    roty(yx_ang, rmtx);
    multi_matrix(orien, 3, 3, rmtx, 3, 3, temp);
    rotz(0.5 * param[6] + phi, rmtx);
    multi_matrix(temp, 3, 3, rmtx, 3, 3, orien);

    /* middle frame orientation w.r.t. 1 */
    rotz(0.5 * param[6] - phi, mst);
    roty(0.5 * yx_ang, rmtx);
    multi_matrix(mst, 3, 3, rmtx, 3, 3, temp);
    rotz(phi, rmtx);
    multi_matrix(temp, 3, 3, rmtx, 3, 3, mst);

    /* base 2 (base-pair 2) position w.r.t. 1 */
    multi_vec_Tmatrix(param, 3, mst, 3, 3, pos);

    /* middle frame position w.r.t. 1 */
    for (i = 1; i <= 3; i++)
        mpos[i] = 0.5 * pos[i];

    free_dmatrix(rmtx, 1, 3, 1, 3);
    free_dmatrix(temp, 1, 3, 1, 3);
}

/* get the orientation and position of one base (base-pair) relative
   to the other given 6 local helical parameters. the middle helical
   frame is also defined
   % ---------------------------------------------------------------
   %   param: 6 base-pair or step parameters in the order of:
   %            x       y       z          x        y       z
   %          x-Disp  y-Disp   Rise   Inclination  Tip    Twist
   % ---------------------------------------------------------------
   %   orien: base 2 (base-pair 2) orientation w.r.t. 1
   %   mst: middle helical frame orientation w.r.t. 1
   %   pos: base 2 (base-pair 2) origin w.r.t. 1
   %   mpos: middle helical frame origin w.r.t. 1
   % --------------------------------------------------------------- */
void xhelfunc(double *param, double **orien, double **mst, double *pos, double *mpos)
{
    double phi, yx_ang;
    double hinge[4];
    double y[4] = { EMPTY_NUMBER, 0.0, 1.0, 0.0 };
    double z[4] = { EMPTY_NUMBER, 0.0, 0.0, 1.0 };
    double htwist;
    double **dxdy, **rmtx, **rot1_h, **rot2_h, **temp;
    long i, j;

    htwist = param[6];

    /* get the yx_ang angle (POSITIVE) and phi angle (-180, +180) */
    hinge[1] = param[4];
    hinge[2] = param[5];
    hinge[3] = 0.0;
    yx_ang = veclen(hinge);
    phi = vec_ang(hinge, y, z);

    /* special case: helical twist = 0.0 */
    if (fabs(htwist) < HTWIST0)  /* XEPS cf: helical_par */
        htwist = HTWIST0;

    rot1_h = dmatrix(1, 3, 1, 3);
    rot2_h = dmatrix(1, 3, 1, 3);
    rmtx = dmatrix(1, 3, 1, 3);
    temp = dmatrix(1, 3, 1, 3);

    /* base 1 (base-pair 1) helical frame w.r.t. 1 */
    rotz(-phi, rot1_h);
    roty(-yx_ang, rmtx);
    multi_matrix(rot1_h, 3, 3, rmtx, 3, 3, temp);
    rotz(phi, rmtx);
    multi_matrix(temp, 3, 3, rmtx, 3, 3, rot1_h);

    /* base 2 (base-pair 2) helical frame w.r.t. 1 */
    rotz(htwist, rmtx);
    multi_matrix(rot1_h, 3, 3, rmtx, 3, 3, rot2_h);

    /* base 2 (base-pair 2) orientation w.r.t. 1 */
    rotz(-phi, orien);
    roty(-yx_ang, rmtx);
    multi_matrix(orien, 3, 3, rmtx, 3, 3, temp);
    rotz(htwist, rmtx);
    multi_matrix(temp, 3, 3, rmtx, 3, 3, orien);
    roty(yx_ang, rmtx);
    multi_matrix(orien, 3, 3, rmtx, 3, 3, temp);
    rotz(phi, rmtx);
    multi_matrix(temp, 3, 3, rmtx, 3, 3, orien);

    /* middle helical frame orientation w.r.t. 1 */
    rotz(0.5 * htwist, rmtx);
    multi_matrix(rot1_h, 3, 3, rmtx, 3, 3, mst);

    dxdy = dmatrix(1, 2, 1, 3);

    /* base 2 (base-pair 2) position w.r.t. 1 */
    for (i = 1; i <= 2; i++)
        for (j = 1; j <= 3; j++)
            dxdy[i][j] = rot2_h[j][i] - rot1_h[j][i];
    multi_vec_matrix(param, 2, dxdy, 2, 3, pos);
    for (i = 1; i <= 3; i++)
        pos[i] += param[3] * rot1_h[i][3];

    /* middle helical frame position w.r.t. 1 */
    for (i = 1; i <= 2; i++)
        for (j = 1; j <= 3; j++)
            dxdy[i][j] = -rot1_h[j][i];
    multi_vec_matrix(param, 2, dxdy, 2, 3, mpos);
    for (i = 1; i <= 3; i++)
        mpos[i] += 0.5 * param[3] * rot1_h[i][3];

    free_dmatrix(rot1_h, 1, 3, 1, 3);
    free_dmatrix(rot2_h, 1, 3, 1, 3);
    free_dmatrix(rmtx, 1, 3, 1, 3);
    free_dmatrix(temp, 1, 3, 1, 3);
    free_dmatrix(dxdy, 1, 2, 1, 3);
}

/* origin xyz coordinates followed by direction cosines of x-, y- & z-axes */
void print_ref_frames(FILE * fp, double *pos_next, double **orien_next)
{
    char *format = "%10.4f%10.4f%10.4f\n";
    long i;

    fprintf(fp, format, pos_next[1], pos_next[2], pos_next[3]);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, orien_next[1][i], orien_next[2][i], orien_next[3][i]);
}
