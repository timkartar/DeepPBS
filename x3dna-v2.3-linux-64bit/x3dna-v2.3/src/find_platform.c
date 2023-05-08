#include "x3dna.h"

typedef struct {
    char pdbfile[BUF512];
    char pair[BUF512];
    char cnt[BUF512];
    char sp_chk[BUF512];  /* check for interactions with a specific pair */

    double overlap;
    double metal;

    long bb_hbfile;
    long gap;
    long min_hb;
    long waters;
} args_pf;

static void set_default_args(args_pf * args)
{
    strcpy(args->pdbfile, "");
    strcpy(args->pair, "");
    strcpy(args->cnt, "ALL");
    strcpy(args->sp_chk, "");

    args->overlap = 0.5;
    args->metal = 0.0;

    args->bb_hbfile = FALSE;
    args->gap = 0;
    args->min_hb = 1;
    args->waters = FALSE;
}

static void platform_usage(void)
{
    if (Gvars.DEBUG > DEBUG_LEVEL)
        fprintf(stderr,
                "Usage: platform [-b] [-gap=INTEGER] [-pair=STRING] [-min_hb=INTEGER] \n"
                "             [-cnt=STRING] [-overlap=FLOAT] [-water] [-metal=cutoff] \n"
                "             [-sp_chk=STRING] pdbfile\n");
    else
        fprintf(stderr, "Usage: find_platform pdbfile\n");
    contact_msg(0);
}

static char *PO[] = { " O1P", " O2P", " O3'", " O4'", " O5'" };

static long numPO = sizeof PO / sizeof PO[0] - 1;

static void get_bb_hbond(long i, long j, long *idx, long **seidx, double **xyz,
                         char **AtomName, char *bb_hbinfo, miscPars * misc_pars)
{
    char *ma, *na;
    long m, n;
    double dval0 = XBIG, dval;

    strcpy(bb_hbinfo, "[**** - **** (*.**)]");  /* default */

    for (m = seidx[i][1]; m <= seidx[i][2]; m++) {
        ma = AtomName[m];
        for (n = seidx[j][1]; n <= seidx[j][2]; n++) {
            na = AtomName[n];
            if (!is_baseatom(ma) && !is_baseatom(na) &&
                (idx[m] == 2 || idx[m] == 4) && (idx[n] == 2 || idx[n] == 4)
                && within_limits(xyz[m], xyz[n], misc_pars->hb_lower, misc_pars->hb_dist1)) {

                if (num_strmatch(ma, PO, 0, numPO) && num_strmatch(na, PO, 0, numPO))
                    continue;  /* no H-bond within PO4 group & O4' */

                dval = p1p2_dist(xyz[m], xyz[n]);
                if (dval < dval0) {
                    dval0 = dval;
                    sprintf(bb_hbinfo, "[%s - %s (%4.2f)]", ma, na, dval);
                }
            }
        }
    }
}

/* '1jbr 01 C:10  C:11'; separated by white space */
static void validate_sp_chk(char *sp_chk, char *pdbfile)
{
    char str0[BUF512], str[BUF512], tmp[BUF512], *items[BUF512];
    long nitem, num_max = 4;
    long LMAX_CUTOFF = 2108888888;  /* 2147483647 */

    if (is_equal_string(sp_chk, ""))
        return;

    upperstr(sp_chk);
    strcpy(str0, sp_chk);  /* do not change str0! */

    nitem = item_list(str0, items, num_max, " \t\n;,");

    if (nitem != num_max)
        fatal("[%s]: wrong format: 4 items as {1jbr 01 C:10  C:11}\n", sp_chk);

    sprintf(str, "%s.PDB", items[1]);
    strcpy(tmp, pdbfile);
    upperstr(tmp);

    if (!is_equal_string(str, tmp))
        fatal("[%s]: wrong format {PDB id as the 1st entry}\n", sp_chk);

    if (cvt2long(items[2]) > LMAX_CUTOFF)
        fatal("[%s]: invalid platform serial number\n", sp_chk);

    if ((items[3][1] != ':') || (items[4][1] != ':'))
        fatal("[%s]: invalid residue identification, e.g., C:10\n", sp_chk);
}

static void validate_pair_str(char *pair)
{
    char *valid_chars = "ACGTUacgtu-+*", *bases = "ACGTUacgtu";

    if (is_empty_string(pair))  /* just '-p=' */
        strcpy(pair, "G-**+-U");

    if (strlen(pair) != 7) {
        fprintf(stderr, "Wrong pair format <%s>: must be 7-char long e.g., 'G-**+-U'\n", pair);
        fatal("Have a look of find_pair/analyze output to see more examples\n");
    }

    if (!string_contains_only_those_characters(pair, valid_chars))
        fatal("With invalid chars <%s> must be [%s]\n", pair, valid_chars);

    if (!strchr(bases, pair[0]) || !strchr(bases, pair[6]))
        fatal("Either 1st or 7th char is not a base <%s> [%s]\n", pair, bases);

    if (pair[1] != '-' || pair[5] != '-')
        fatal("The 2nd/6th char must be '-' <%s>\n", pair);

    if (pair[2] != '-' && pair[2] != '*')
        fatal("The 3rd char must be [-*] for WC/non-WC pair <%s>\n", pair);

    if (pair[3] != '-' && pair[3] != '*')
        fatal("The 4th char must be [-*] for WC/non-WC geometry <%s>\n", pair);

    if (pair[4] != '-' && pair[4] != '+')
        fatal("The 5th char must be [-+] for anti/parallel z-axes <%s>\n", pair);
}

static void platform_cmdline(int argc, char *argv[], args_pf * args)
{
    long i;

    if (argc < 2)
        platform_usage();

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        if (case_str_pmatch(argv[i], "-b"))
            args->bb_hbfile = TRUE;

        else if (case_str_pmatch(argv[i], "-g"))
            args->gap = get_lvalue(argv[i], 0, 5);

        else if (case_str_pmatch(argv[i], "-mi"))
            args->min_hb = get_lvalue(argv[i], 1, BUF512);

        else if (case_str_pmatch(argv[i], "-w"))
            args->waters = TRUE;

        else if (case_str_pmatch(argv[i], "-me")) {
            if (!strchr(argv[i], '='))  /* just -me: use default cut-off */
                args->metal = 6.0;
            else
                args->metal = get_dvalue(argv[i], 0.0, 20.0);

        } else if (case_str_pmatch(argv[i], "-s"))
            get_strvalue(argv[i], args->sp_chk, FALSE);

        else if (case_str_pmatch(argv[i], "-c"))
            get_strvalue(argv[i], args->cnt, FALSE);

        else if (case_str_pmatch(argv[i], "-p")) {
            if (!strchr(argv[i], '='))  /* just -p: use default G+U platform pair */
                strcpy(args->pair, "G-**+-U");
            else {
                get_strvalue(argv[i], args->pair, FALSE);
                validate_pair_str(args->pair);
            }

        } else if (case_str_pmatch(argv[i], "-o"))
            args->overlap = get_dvalue(argv[i], 0, 5);

        else
            platform_usage();
    }

    if (argc == i + 1) {
        strcpy(args->pdbfile, argv[i]);
        if (!exist_file(args->pdbfile))
            fatal("PDB file <%s> does not exist\n", args->pdbfile);
    } else
        platform_usage();

    upperstr(args->cnt);  /* case does not matter here */

    validate_sp_chk(args->sp_chk, args->pdbfile);

    if ((is_empty_string(args->pair) && is_empty_string(args->sp_chk)) && args->metal) {
        fprintf(stderr, "-metal option must be with -pair: ignore %g\n", args->metal);
        args->metal = 0.0;
    }
}

static void match_pair(long min_hb, char *pdbid, char *pair, char *parfile)
{
    char *p0, *p1, *line, par[BUF512], hb[BUF512], pair_id[BUF512];
    long iser, i, j, num_hb, num_mp = 0;
    FILE *fp_par, *fp_o2s;

    remove_file(O2_STACK);
    if (is_empty_string(pair) || !exist_file(parfile))
        return;

    fp_o2s = open_file(O2_STACK, "w");

    fprintf(fp_o2s, "# <id>%s</id>\n", pdbid);
    fprintf(fp_o2s, "# <pair>%s</pair>\n", pair);

    fp_par = open_file(parfile, "r");

    while ((p0 = my_getline(fp_par)) != NULL) {
        line = trim(p0);  /* keep the original value of p0 */
        if (line[0] == '#' || line[0] == '\0' || case_str_pmatch(line, "Base-pair criteria used")) {
            free(p0);
            continue;
        }

        /* match line "1    9   10 ....>A:2655_:[..G]G-**+-U[..U]:2656_:A<...." */
        if (sscanf(line, "%ld %ld %ld %s", &iser, &i, &j, pair_id) == 4) {
            if (strchr(line, ']') == NULL)
                fatal("wrong format for pair-identity <%s>\n", line);

            p1 = my_getline(fp_par);  /* G+U 7.77    1.83   -0.27    0.82    5.94   12.44 */
            line = trim(p1);
            strcpy(par, line + 5);
            free(p1);

            p1 = my_getline(fp_par);  /* G+U [3]  O2'- O1P 3.16  O2'- O2P 2.73  N2 - O4  2.64 */
            line = trim(p1);
            strcpy(hb, line);
            free(p1);

            if (sscanf(hb, "%*s [%ld]", &num_hb) != 1)
                fatal("wrong format for H-bonding info <%s>\n", hb);

            if (case_str_pmatch(pair_id + 18, pair) && num_hb >= min_hb) {
                fprintf(stderr, "num_hb: %ld\n", num_hb);
                num_mp++;
                fprintf(fp_o2s, "%5ld %5ld %5ld %s %s {%s}\n", num_mp, i, j, pair_id, hb, par);
            }
        }

        free(p0);
    }

    close_file(fp_o2s);
    close_file(fp_par);

    if (!num_mp)
        remove_file(O2_STACK);
}

static void get_gap_name(char *bname, long gap, char *ext, char *ename)
{
    if (gap)
        sprintf(ename, "%s_%ld.%s", bname, gap, ext);
    else  /* gap=0 is the default dinucleotide platform */
        sprintf(ename, "%s.%s", bname, ext);
}

static void P_C4_parameters(long ia, long ib, long **seidx, char **AtomName, double **xyz,
                            FILE * fp_par)
{
    long i, k = 0;
    double d[5], v1[4], v2[4], **xyz4;  /* P(ia), C4'(ia), P(ib), C4'(ib) */

    xyz4 = dmatrix(1, 4, 1, 3);

    for (i = seidx[ia][1]; i <= seidx[ia][2]; i++) {
        if (is_equal_string(AtomName[i], " P  ")) {
            copy_dvector(xyz4[1], xyz[i], 1, 3);
            k++;
        } else if (is_equal_string(AtomName[i], " C4'")) {
            copy_dvector(xyz4[2], xyz[i], 1, 3);
            k++;
        }
    }

    for (i = seidx[ib][1]; i <= seidx[ib][2]; i++) {
        if (is_equal_string(AtomName[i], " P  ")) {
            copy_dvector(xyz4[3], xyz[i], 1, 3);
            k++;
        } else if (is_equal_string(AtomName[i], " C4'")) {
            copy_dvector(xyz4[4], xyz[i], 1, 3);
            k++;
        }
    }

    if (k != 4)
        fprintf(fp_par, "        NA missing atom P(ia), C4'(ia), P(ib) or C4'(ib)\n");
    else {
        d[0] = torsion2(xyz4);  /* torsion P(ia), C4'(ia), P(ib), C4'(ib) */

        ddxyz(xyz4[2], xyz4[1], v1);  /* C4'(ia)-->P(ia) */
        ddxyz(xyz4[2], xyz4[3], v2);  /* C4'(ia)-->P(ib) */
        d[1] = magang(v1, v2);

        ddxyz(xyz4[3], xyz4[2], v1);  /* P(ib)-->C4'(ia) */
        ddxyz(xyz4[3], xyz4[4], v2);  /* P(ib)-->C4'(ib) */
        d[2] = magang(v1, v2);

        d[3] = p1p2_dist(xyz4[1], xyz4[3]);  /* P(ia)--P(ib) */
        d[4] = p1p2_dist(xyz4[2], xyz4[4]);  /* C4'(ia)--C4'(ib) */

        fprintf(fp_par, "        P(i)-C4'(i)-P(i+1)-C4'(i+1): "
                "%8.0f %8.0f %8.0f %8.1f %8.1f\n", d[0], d[1], d[2], d[3], d[4]);
    }

    free_dmatrix(xyz4, 1, 4, 1, 3);
}

static void check_neighbor_pair(long num, long num_residue, long *RY, long *idx,
                                double **NC1xyz, double **o3_p, double **orien,
                                double **org, miscPars * misc_pars, long **seidx,
                                double **xyz, char **AtomName, char **ResName,
                                char *ChainID, long *ResSeq, char **Miscs, char *bseq,
                                long gap, long min_hb, long waters, char *pdbid, char *pair)
{
    char wc[BUF512], pair_id[BUF512], b1[BUF512], b2[BUF512], str[BUF512];
    char strfile[BUF512], parfile[BUF512];
    double morg[4], morien[10];
    double rtn_val[RTNNUM], *chi;
    long bpid, i, ir, j, jr, **ring_atom, **htm_water;
    long inum_base = 2, num_bp = 0, ivec[3];
    FILE *fp_par, *fp_str;

    chi = dvector(1, num_residue);
    get_chi_angle(num_residue, RY, bseq, seidx, xyz, AtomName, ResName, ChainID,
                  ResSeq, Miscs, chi, NULL);

    htm_water = lmatrix(1, 4, 0, num);  /* HETATM and water index */
    init_htm_water(waters, num, num_residue, idx, htm_water);
    identify_htw(num_residue, seidx, RY, AtomName, ResName, ChainID, ResSeq,
                 Miscs, xyz, htm_water);

    ring_atom = lmatrix(1, num_residue, 1, 19);
    ring_oidx(num, num_residue, RY, seidx, AtomName, xyz, idx, ring_atom);

    get_gap_name(pdbid, gap, "str", strfile);
    fp_str = open_file(strfile, "w");

    get_gap_name(pdbid, gap, "par", parfile);
    fp_par = open_file(parfile, "w");

    print_bp_crit(misc_pars, fp_par);

    for (i = 1; i < num_residue - gap; i++) {
        if (RY[i] < 0)
            continue;

        j = i + 1 + gap;  /* gap = 0 for i --> i + 1 */
        if (RY[j] < 0 || !is_linked_by_gap(i, j, o3_p))
            continue;

        check_pair(i, j, bseq, seidx, xyz, NC1xyz, orien, org, idx, AtomName,
                   misc_pars, rtn_val, &bpid, ring_atom, 0);

        if (!bpid)
            continue;

        num_bp++;
        bpid_wc_str(bpid, rtn_val[35], wc);
        sprintf(pair_id, "%c-%s-%c", bseq[i], wc, bseq[j]);

        ir = seidx[i][1];
        jr = seidx[j][1];
        base_str(ChainID[ir], ResSeq[ir], Miscs[ir], ResName[ir], bseq[i], 1, b1);
        base_str(ChainID[jr], ResSeq[jr], Miscs[jr], ResName[jr], bseq[j], 2, b2);

        get_bb_hbond(i, j, idx, seidx, xyz, AtomName, str, misc_pars);
        fprintf(fp_par, "%5ld %5ld %5ld %s-%s-%s %s\n", num_bp, i, j, b1, wc, b2, str);

        print_pairinfo(i, j, bseq[i], bseq[j], rtn_val, chi, misc_pars, seidx, idx,
                       AtomName, xyz, bseq, 1, fp_par);

        P_C4_parameters(i, j, seidx, AtomName, xyz, fp_par);

        /* multiple structure file containing all base-pairs */
        sprintf(str, "%s-%s-%s", b1, wc, b2);
        fprintf(fp_str, "%6s    %4ld\n", "MODEL ", num_bp);
        fprintf(fp_str, "REMARK    Section #%4.4ld %s\n", num_bp, str);
        fprintf(fp_str, "REMARK    %s\n", Gvars.X3DNA_VER);
        ivec[1] = i;
        ivec[2] = j;
        pair2mst(inum_base, ivec, AtomName, ResName, ChainID, ResSeq, Miscs, xyz,
                 orien, org, seidx, morien, morg, htm_water, misc_pars, fp_str);
        fprintf(fp_str, "ENDMDL\n");
    }

    free_dvector(chi, 1, num_residue);
    free_lmatrix(ring_atom, 1, num_residue, 1, 19);
    free_lmatrix(htm_water, 1, 4, 0, num);

    close_file(fp_par);
    close_file(fp_str);

    if (!num_bp) {
        remove_file(strfile);
        remove_file(parfile);
    }

    match_pair(min_hb, pdbid, pair, parfile);
}

/* calculate all possible H-bonds between neighbouring residues: an
 * independent function() with option -b for 'bb_hbfile' */
static void check_neighbor_hbond(long num_residue, long *RY, long *idx, double **o3_p,
                                 long **seidx, double **xyz, char **AtomName,
                                 char **ResName, char *ChainID, long *ResSeq,
                                 char **Miscs, char *bseq, char *pdbid, long bb_hbfile,
                                 long gap, miscPars * misc_pars)
{
    char *ma, *na, b1[BUF512], b2[BUF512], hbdfile[BUF512];
    long i, j, m, n, num_hb = 0;
    double dval;

    FILE *fp;

    if (!bb_hbfile)
        return;

    get_gap_name(pdbid, gap, "hbd", hbdfile);
    fp = open_file(hbdfile, "w");

    for (i = 1; i < num_residue - gap; i++) {
        if (RY[i] < 0)
            continue;

        j = i + 1 + gap;  /* gap = 0 for i --> i + 1 */
        if (RY[j] < 0 || !is_linked_by_gap(i, j, o3_p))
            continue;

        for (m = seidx[i][1]; m <= seidx[i][2]; m++) {
            ma = AtomName[m];
            for (n = seidx[j][1]; n <= seidx[j][2]; n++) {
                na = AtomName[n];
                if (!is_baseatom(ma) && !is_baseatom(na) &&
                    (idx[m] == 2 || idx[m] == 4) && (idx[n] == 2 || idx[n] == 4)
                    && within_limits(xyz[m], xyz[n], misc_pars->hb_lower, misc_pars->hb_dist1)) {

                    if (num_strmatch(ma, PO, 0, numPO) && num_strmatch(na, PO, 0, numPO))
                        continue;  /* no H-bond within PO4 group & O4' */

                    num_hb++;

                    base_str(ChainID[m], ResSeq[m], Miscs[m], ResName[m], bseq[i], 1, b1);
                    base_str(ChainID[n], ResSeq[n], Miscs[n], ResName[n], bseq[j], 1, b2);

                    dval = p1p2_dist(xyz[m], xyz[n]);
                    fprintf(fp, "%s\t%ld\t%ld\t%s\t%s", pdbid, i, j, b1, b2);
                    fprintf(fp, "\t%s\t%s\t%4.2f\t%ld\t%ld\n", ma, na, dval, m, n);
                }
            }
        }
    }

    close_file(fp);

    if (!num_hb)
        remove_file(hbdfile);
}

static long match_residue_serial_number(char *res, long num_residue, long **seidx,
                                        char **ResName, char *bseq, char *ChainID, long *ResSeq)
{
    char str[BUF512];
    long i, k, num = 0, idx = 0;

    for (i = 1; i <= num_residue; i++) {
        k = seidx[i][1];
        if (is_equal_string(ResName[k], "HOH"))
            continue;  /* special case (3cc2): HETATM96435  O   HOH 01370: G 01370 */

        sprintf(str, "%c:%c%ld", ChainID[k], bseq[i], ResSeq[k]);
        if (is_equal_string(res, str)) {
            num++;
            idx = i;
        }
    }

    if (num == 0)
        fatal("no matches for residue: %s\n", res);
    else if (num > 1)
        fatal("%ld matches for residue: %s\n", num, res);

    return idx;
}

static void check_hb_connection(long ir, long num_residue, long **seidx, double **xyz,
                                char **AtomName, char **ResName, long *idx,
                                miscPars * misc_pars, double **orien, double **org,
                                long *inum, long *ivec)
{
    long i, ib, ie, j, k, m = *inum, b1, b2;
    double dv, dd[4], HB_UPPER = 3.2, HDV = 2.0;

    for (i = seidx[ir][1]; i <= seidx[ir][2]; i++) {
        if (idx[i] != 2 && idx[i] != 4)  /* non N/O */
            continue;

        b1 = is_baseatom(AtomName[i]);

        for (j = 1; j <= num_residue; j++) {
            if (lval_in_set(j, 1, m, ivec))  /* already counted */
                continue;

            ib = seidx[j][1];
            ie = seidx[j][2];
            if (ie - ib + 1 > 1 && is_equal_string("HOH", ResName[ib]))
                continue;  /* 3cc2: HOH 01058 etc has many waters */

            for (k = seidx[j][1]; k <= seidx[j][2]; k++) {
                if (idx[k] != 2 && idx[k] != 4)  /* non N/O */
                    continue;

                if (num_strmatch(AtomName[i], PO, 0, numPO) &&
                    num_strmatch(AtomName[k], PO, 0, numPO))
                    continue;  /* no H-bond within PO4 group & O4' */

                if (!within_limits(xyz[i], xyz[k], misc_pars->hb_lower, HB_UPPER))
                    continue;

                b2 = is_baseatom(AtomName[k]);

                if (!b1 && !b2) {  /* within backbone atoms */
                    ivec[++m] = j;
                    break;
                }

                if (b1) {  /* with base atom on ir */
                    ddxyz(xyz[k], org[ir], dd);
                    dv = fabs(dot(dd, &orien[ir][6]));
                    if (dv <= HDV) {
                        ivec[++m] = j;
                        break;
                    }
                }

                if (b2) {  /* with base atom on j */
                    ddxyz(xyz[i], org[j], dd);
                    dv = fabs(dot(dd, &orien[j][6]));
                    if (dv <= HDV) {
                        ivec[++m] = j;
                        break;
                    }
                }
            }
        }
    }

    *inum = m;
}

/* only with ivec[1] and [2], i.e., the first two residue in platform pair */
static void check_metal(double cutoff, long num_residue, long **seidx, double **xyz,
                        char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                        long *idx, long *inum, long *ivec, FILE * fp)
{
    char b1[BUF512], b2[BUF512];
    double dval, mdist;
    long i, im, ires, j, k, nchk = 2, m = *inum;
    long midx, num_atoms, *is_metal;

    if (cutoff == 0.0)
        return;

    num_atoms = seidx[num_residue][2];
    is_metal = lvector(1, num_atoms);
    atom_metal(num_atoms, AtomName, is_metal);

    for (i = 1; i <= num_residue; i++) {
        if (lval_in_set(i, 1, m, ivec))  /* already counted */
            continue;

        for (im = seidx[i][1]; im <= seidx[i][2]; im++) {
            if (!is_metal[im])  /* skip non-metal atoms */
                continue;

            midx = FALSE;
            mdist = XBIG;

            for (j = 1; j <= nchk; j++) {
                ires = ivec[j];

                for (k = seidx[ires][1]; k <= seidx[ires][2]; k++) {
                    if (idx[k] != 2 && idx[k] != 4)  /* non N/O */
                        continue;

                    dval = p1p2_dist(xyz[im], xyz[k]);
                    if (dval <= cutoff && dval < mdist) {
                        mdist = dval;
                        midx = k;
                    }
                }
            }

            if (midx) {
                ivec[++m] = i;

                residue_idstr(ChainID[im], ResSeq[im], ResName[im], b1);
                residue_idstr(ChainID[midx], ResSeq[midx], ResName[midx], b2);

                fprintf(stderr, "\tMetal %s [%s] with residue %s [%s] (%.2f)\n",
                        b1, AtomName[im], b2, AtomName[midx], mdist);
                fprintf(fp, "REMARK    Metal: %s [%s] with residue %s [%s] (%.2f)\n",
                        b1, AtomName[im], b2, AtomName[midx], mdist);
            }
        }
    }

    free_lvector(is_metal, 1, num_atoms);

    *inum = m;
}

static void check_stacking(long ir, long num_residue, double **xyz, double **orien,
                           double **org, long **ring_atom, miscPars * misc_pars,
                           long *inum, long *ivec, FILE * fp, double overlap)
{
    long i, m = *inum;
    double dv, dd[4], oave[4], zave[4], OV = overlap;

    for (i = 1; i <= num_residue; i++) {
        if (lval_in_set(i, 1, m, ivec))  /* already counted */
            continue;

        get_bp_zoave(ir, i, orien, org, oave, zave);

        ddxyz(org[ir], org[i], dd);
        dv = fabs(dot(dd, zave));

        if ((veclen(dd) > misc_pars->max_dorg) || !dval_in_range(dv, 2.0, 5.0))
            continue;

        dd[1] = get_oarea(ir, i, ring_atom, oave, zave, xyz, 0);
        dd[2] = get_oarea(ir, i, ring_atom, org[ir], &orien[ir][6], xyz, 0);
        dd[3] = get_oarea(ir, i, ring_atom, org[i], &orien[i][6], xyz, 0);

        if (dd[1] > OV || dd[2] > OV || dd[3] > OV) {
            m++;
            ivec[m] = i;
            fprintf(stderr, "\t%ld vs %ld\t%5.2f\t%5.2f\t%5.2f\n", ir, i, dd[1], dd[2], dd[3]);
            fprintf(fp, "REMARK    Overlap: %4ld vs %4ld  m=%5.2f  i=%5.2f  j=%5.2f\n",
                    ir, i, dd[1], dd[2], dd[3]);
        }
    }

    *inum = m;
}

/* check specifically for interactions with O6 and N6 of G in G+U platform */
static void check_gu_hoogsteen(char *bp, long ir, long num_residue, long **seidx,
                               double **xyz, char **AtomName, char **ResName,
                               char *ChainID, long *ResSeq, long *idx,
                               miscPars * misc_pars, double **orien, double **org,
                               long inum, long *ivec, FILE * fp)
{
    char b1[BUF512], b2[BUF512];
    long i, j, k = 0, midx;
    double dval, dv, mdist, dd[4], HB_UPPER = 3.2, HDV = 2.0;

    if (!is_equal_string(bp, "G-U"))
        return;

    for (i = seidx[ir][1]; i <= seidx[ir][2]; i++) {
        if (!is_equal_string(AtomName[i], " O6 ")
            && !is_equal_string(AtomName[i], " N7 "))
            continue;

        midx = FALSE;
        mdist = XBIG;

        for (j = 1; j <= num_residue; j++) {
            if (!lval_in_set(j, 3, inum, ivec))  /* should be already counted */
                continue;

            for (k = seidx[j][1]; k <= seidx[j][2]; k++) {
                if (idx[k] != 2 && idx[k] != 4)  /* non N/O */
                    continue;

                ddxyz(xyz[k], org[ir], dd);
                dv = fabs(dot(dd, &orien[ir][6]));
                if (dv > HDV)
                    continue;

                dval = p1p2_dist(xyz[i], xyz[k]);
                if (!dval_in_range(dval, misc_pars->hb_lower, HB_UPPER))
                    continue;

                if (dval < mdist) {
                    mdist = dval;
                    midx = k;
                }
            }
        }

        if (midx) {
            residue_idstr(ChainID[i], ResSeq[i], ResName[i], b1);
            residue_idstr(ChainID[midx], ResSeq[midx], ResName[midx], b2);

            fprintf(stderr, "\t%s [%s] with residue %s [%s] (%.2f)\n", b1, AtomName[i],
                    b2, AtomName[midx], mdist);
            fprintf(fp, "REMARK    %s [%s] with residue %s [%s] (%.2f)\n", b1,
                    AtomName[i], b2, AtomName[midx], mdist);
        }
    }
}

static void output_residue_list(long inum, long *ivec, long **seidx, char **ResName,
                                char *bseq, char *ChainID, long *ResSeq, FILE * fp)
{
    char cid, cnt[BUF512] = "", str[BUF512], res[BUF512];
    long i, j, k, n = 0, ib;

    for (i = 1; i <= inum; i++) {
        k = ivec[i];
        ib = seidx[k][1];
        if (!is_equal_string(ResName[ib], "HOH")) {  /* skip water */
            n++;
            fprintf(stderr, "%4ld\t%ld\t%c\t%4ld\t%s\n", n, k, ChainID[ib],
                    ResSeq[ib], ResName[ib]);
            cid = ChainID[ib];
            if (cid == ' ')
                cid = '_';
            if (bseq[k] == '\0') {  /* for non-nucleotide */
                sprintf(res, "%s", ResName[ib]);
                for (j = 0; j <= 2; j++)
                    if (res[j] == ' ')
                        res[j] = '_';
            } else
                sprintf(res, "%c", bseq[k]);
            sprintf(str, "%c:%s%ld  ", cid, res, ResSeq[ib]);
            strcat(cnt, str);
        }
    }

    fprintf(fp, "REMARK    Contacts: %s\n", cnt);
}

/* write a set of 'inum' residues in ivev[] with reference to the middle-frame
 * of the pair defined by the first two, ivec[1], ivec[2] */
static void reorient_pair(long inum, long *ivec, char **AtomName, char **ResName,
                          char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                          double **orien, double **org, long **seidx, FILE * fp)
{
    double morg[4], **mst, **xyz_residue;
    long i, ik, j, m, iser = 0;

    mst = dmatrix(1, 3, 1, 3);
    xyz_residue = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);

    init_dvector(morg, 1, 3, 0.0);

    cehs_average(2, ivec, orien, org, mst, morg);

    for (i = 1; i <= inum; i++) {
        ik = ivec[i];
        for (j = seidx[ik][1]; j <= seidx[ik][2]; j++) {
            m = j - seidx[ik][1] + 1;
            cpxyz(xyz[j], xyz_residue[m]);
        }
        change_xyz(0, morg, mst, seidx[ik][2] - seidx[ik][1] + 1, xyz_residue);
        pdb_record(seidx[ik][1], seidx[ik][2], &iser, 1, AtomName, ResName,
                   ChainID, ResSeq, xyz_residue, Miscs, fp);
    }

    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(xyz_residue, 1, NUM_RESIDUE_ATOMS, 1, 3);
}

/* private and specific: check for tertiary interactions with a
 * specific pair: e.g., a G+U with 2 H-bonds */
static void check_specific_pair(char *pdbid, long num, long num_residue, double **orien,
                                double **org, miscPars * misc_pars, long **seidx,
                                long *idx, double **xyz, char **AtomName, char **ResName,
                                char *ChainID, long *ResSeq, char **Miscs, long *RY,
                                char *bseq, char *sp_chk, double overlap, double metal)
{
    char res1[BUF512], res2[BUF512], cagfile[BUF512], *items[BUF512];
    char b1[BUF512], b2[BUF512], bp[BUF512];
    long i, iser, j, k, inum, num_max = 4, ivec[BUF512], **ring_atom;
    FILE *fp;

    if (is_equal_string(sp_chk, ""))
        return;

    item_list(sp_chk, items, num_max, " \t\n;,");
    iser = cvt2long(items[2]);
    strcpy(res1, items[3]);
    strcpy(res2, items[4]);

    i = match_residue_serial_number(res1, num_residue, seidx, ResName, bseq, ChainID, ResSeq);
    j = match_residue_serial_number(res2, num_residue, seidx, ResName, bseq, ChainID, ResSeq);
    sprintf(bp, "%c-%c", bseq[i], bseq[j]);  /* e.g., G-U */

    inum = 2;
    ivec[1] = i;
    ivec[2] = j;

    k = seidx[i][1];
    base_str(ChainID[k], ResSeq[k], Miscs[k], ResName[k], bseq[i], 1, b1);
    k = seidx[j][1];
    base_str(ChainID[k], ResSeq[k], Miscs[k], ResName[k], bseq[j], 1, b2);

    sprintf(cagfile, "%s_%2.2ld.cag", pdbid, iser);
    fp = open_file(cagfile, "w");

    fprintf(fp, "REMARK    PDB_ID: %s\n", pdbid);
    fprintf(fp, "REMARK    RES1: %s [%s]\n", res1, b1);
    fprintf(fp, "REMARK    RES2: %s [%s]\n", res2, b2);
    fprintf(fp, "REMARK    ISER: %ld\n", iser);
    fprintf(fp, "REMARK    ");
    print_sep(fp, '-', 68);

    ring_atom = lmatrix(1, num_residue, 1, 19);
    ring_oidx(num, num_residue, RY, seidx, AtomName, xyz, idx, ring_atom);

    fprintf(fp, "REMARK    Section #%4.4ld\n", iser);
    fprintf(fp, "REMARK    %s\n", Gvars.X3DNA_VER);

    check_stacking(i, num_residue, xyz, orien, org, ring_atom, misc_pars, &inum,
                   ivec, fp, overlap);
    check_stacking(j, num_residue, xyz, orien, org, ring_atom, misc_pars, &inum,
                   ivec, fp, overlap);

    check_hb_connection(i, num_residue, seidx, xyz, AtomName, ResName, idx,
                        misc_pars, orien, org, &inum, ivec);
    check_hb_connection(j, num_residue, seidx, xyz, AtomName, ResName, idx,
                        misc_pars, orien, org, &inum, ivec);

    check_gu_hoogsteen(bp, i, num_residue, seidx, xyz, AtomName, ResName, ChainID,
                       ResSeq, idx, misc_pars, orien, org, inum, ivec, fp);

    check_metal(metal, num_residue, seidx, xyz, AtomName, ResName, ChainID, ResSeq,
                idx, &inum, ivec, fp);

    output_residue_list(inum, ivec, seidx, ResName, bseq, ChainID, ResSeq, fp);

    reorient_pair(inum, ivec, AtomName, ResName, ChainID, ResSeq, Miscs, xyz,
                  orien, org, seidx, fp);

    fprintf(fp, "ENDMDL\n");

    free_lmatrix(ring_atom, 1, num_residue, 1, 19);

    close_file(fp);
}

/* cnt: "ALL", "IALL", "ISTACK", "IHBOND", "JALL", "JSTACK", "JHBOND" */
static void check_contact(char *pdbid, char *pair, long num, long num_residue,
                          double **orien, double **org, miscPars * misc_pars,
                          long **seidx, long *idx, double **xyz, char **AtomName,
                          char **ResName, char *bseq, char *ChainID, long *ResSeq,
                          char **Miscs, long *RY, char *cnt, double overlap, double metal)
{
    char *p0, *line, str[BUF512], cntfile[BUF512];
    long i, j, k, inum, ivec[BUF512], **ring_atom;
    FILE *fp, *fp_cnt;

    if (!exist_file(O2_STACK))
        return;

    get_gap_name(pdbid, 0, "cnt", cntfile);
    fp_cnt = open_file(cntfile, "w");

    fprintf(fp_cnt, "REMARK    PDB_ID: %s\n", pdbid);
    fprintf(fp_cnt, "REMARK    PAIR: %s\n", pair);
    fprintf(fp_cnt, "REMARK    CONTACT: %s\n", cnt);
    fprintf(fp_cnt, "REMARK    ");
    print_sep(fp_cnt, '-', 68);

    ring_atom = lmatrix(1, num_residue, 1, 19);
    ring_oidx(num, num_residue, RY, seidx, AtomName, xyz, idx, ring_atom);

    fp = open_file(O2_STACK, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        line = trim(p0);  /* keep the original value of p0 */
        if (line[0] == '#' || line[0] == '\0') {
            free(p0);
            continue;  /* skip empty and commented lines */
        }

        if (sscanf(line, "%ld %ld %ld %s", &k, &i, &j, str) != 4)
            fatal("wrong format <%s>\n", line);

        inum = 2;
        ivec[1] = i;
        ivec[2] = j;

        print_sep(stderr, '-', 76);

        /* multiple structure file containing all base-pairs */
        fprintf(fp_cnt, "%6s    %4ld\n", "MODEL ", k);
        fprintf(fp_cnt, "REMARK    Section #%4.4ld %s\n", k, str);
        fprintf(fp_cnt, "REMARK    %s\n", Gvars.X3DNA_VER);

        if (str_pmatch(cnt, "A") || str_pmatch(cnt, "IA") || str_pmatch(cnt, "IS"))
            check_stacking(i, num_residue, xyz, orien, org, ring_atom, misc_pars,
                           &inum, ivec, fp_cnt, overlap);

        if (str_pmatch(cnt, "A") || str_pmatch(cnt, "JA") || str_pmatch(cnt, "JS"))
            check_stacking(j, num_residue, xyz, orien, org, ring_atom, misc_pars,
                           &inum, ivec, fp_cnt, overlap);

        if (str_pmatch(cnt, "A") || str_pmatch(cnt, "IA") || str_pmatch(cnt, "IH"))
            check_hb_connection(i, num_residue, seidx, xyz, AtomName, ResName, idx,
                                misc_pars, orien, org, &inum, ivec);

        if (str_pmatch(cnt, "A") || str_pmatch(cnt, "JA") || str_pmatch(cnt, "JH"))
            check_hb_connection(j, num_residue, seidx, xyz, AtomName, ResName, idx,
                                misc_pars, orien, org, &inum, ivec);

        check_metal(metal, num_residue, seidx, xyz, AtomName, ResName, ChainID,
                    ResSeq, idx, &inum, ivec, fp_cnt);

        output_residue_list(inum, ivec, seidx, ResName, bseq, ChainID, ResSeq, fp_cnt);

        reorient_pair(inum, ivec, AtomName, ResName, ChainID, ResSeq, Miscs, xyz,
                      orien, org, seidx, fp_cnt);

        fprintf(fp_cnt, "ENDMDL\n");

        free(p0);
    }

    free_lmatrix(ring_atom, 1, num_residue, 1, 19);

    close_file(fp);
    close_file(fp_cnt);
}

static void platform_pars(args_pf * args)
{
    char pdbid[BUF512], *ChainID, *bseq, **AtomName, **ResName, **Miscs;
    double **xyz, **orien, **org, **NC1xyz, **o3_p;
    long num, hetatm = 1, num_residue;
    long *ResSeq, *RY, *idx, **seidx;

    del_extension(args->pdbfile, pdbid);

    /* read in the PDB file */
    num = number_of_atoms(args->pdbfile, hetatm, Gvars.misc_pars.alt_list);
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(args->pdbfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs,
             hetatm, Gvars.misc_pars.alt_list);

    /* get the numbering information of each residue */
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

    /* get base sequence, RY identification */
    bseq = cvector(1, num_residue);
    RY = lvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);

    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 3);
    NC1xyz = dmatrix(1, num_residue, 1, 7);  /* RN9/YN1 & C1' atomic coordinates */
    o3_p = dmatrix(1, num_residue, 1, 8);  /* O3'/P atomic coordinates */

    base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
              ResSeq, Miscs, xyz, orien, org, NC1xyz, o3_p);

    idx = lvector(1, num);
    atom_idx(num, AtomName, NULL, idx);

    check_neighbor_pair(num, num_residue, RY, idx, NC1xyz, o3_p, orien, org,
                        &Gvars.misc_pars, seidx, xyz, AtomName, ResName, ChainID, ResSeq,
                        Miscs, bseq, args->gap, args->min_hb, args->waters, pdbid, args->pair);

    check_specific_pair(pdbid, num, num_residue, orien, org, &Gvars.misc_pars, seidx, idx,
                        xyz, AtomName, ResName, ChainID, ResSeq, Miscs, RY, bseq,
                        args->sp_chk, args->overlap, args->metal);

    check_contact(pdbid, args->pair, num, num_residue, orien, org, &Gvars.misc_pars,
                  seidx, idx, xyz, AtomName, ResName, bseq, ChainID, ResSeq, Miscs, RY,
                  args->cnt, args->overlap, args->metal);

    check_neighbor_hbond(num_residue, RY, idx, o3_p, seidx, xyz, AtomName, ResName,
                         ChainID, ResSeq, Miscs, bseq, pdbid, args->bb_hbfile, args->gap,
                         &Gvars.misc_pars);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cvector(bseq, 1, num_residue);
    free_lvector(RY, 1, num_residue);
    free_lvector(idx, 1, num);

    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 3);
    free_dmatrix(NC1xyz, 1, num_residue, 1, 7);
    free_dmatrix(o3_p, 1, num_residue, 1, 8);
}

int main(int argc, char *argv[])
{
    args_pf args;
    time_t time0;

    time(&time0);

    set_my_globals(argv[0]);

    set_default_args(&args);
    platform_cmdline(argc, argv, &args);

    fprintf(stderr, "\nhandling file <%s>\n", args.pdbfile);
    platform_pars(&args);

    clear_my_globals();

    print_used_time(time0);

    return 0;
}
