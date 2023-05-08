#include "x3dna.h"

static void rebuild_usage(void)
{
    help3dna_usage("rebuild");
}

static void rebuild_cmdline(int argc, char *argv[], char *options[], char *user_input,
                            char *inpfile, char *outfile, long *neg_xdir, long *xml)
{
    long i, j, nop = 4;

    if (argc < 2)
        rebuild_usage();

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        upperstr(argv[i]);

        *xml = get_xmlArgNumber(argv[i], "-X");
        if (*xml)
            continue;

        if (!strcmp(argv[i], "-NEGX"))
            *neg_xdir = 1;
        else {
            for (j = 0; j < nop; j++)
                if (!strcmp(argv[i], options[j])) {
                    strcpy(user_input, argv[i]);
                    break;
                }
            if (j >= nop)
                rebuild_usage();
        }
    }

    if (argc == i + 1) {
        strcpy(inpfile, argv[i]);
        strcpy(outfile, "stdout");
    } else if (argc == i + 2) {
        strcpy(inpfile, argv[i]);
        strcpy(outfile, argv[i + 1]);
    } else
        rebuild_usage();
}

typedef struct {
    long num_bp;
    long is_single;
    long is_helical;
    long xdir;
    long parallel;
    char **bp_seq;
    double **bp_par;
    double **step_par;
    long *pidx;
} struct_pars;

static void init_pars(struct_pars * pars)
{
    pars->num_bp = FALSE;
    pars->is_single = FALSE;
    pars->is_helical = FALSE;
    pars->xdir = FALSE;
    pars->parallel = FALSE;

    pars->bp_seq = NULL;
    pars->bp_par = NULL;
    pars->step_par = NULL;
    pars->pidx = NULL;
}

static void free_pars(struct_pars * pars)
{
    if (pars->bp_seq)
        free_cmatrix(pars->bp_seq, 1, DUMMY, 0, DUMMY);

    if (pars->bp_par)
        free_dmatrix(pars->bp_par, 1, DUMMY, 1, DUMMY);

    if (pars->step_par)
        free_dmatrix(pars->step_par, 1, DUMMY, 1, DUMMY);

    if (pars->pidx)
        free_lvector(pars->pidx, 1, DUMMY);
}

static void read_pars(char *inpfile, struct_pars * pars)
{
    char str[BUF512];
    char *p, *item[14];
    double twist = 0, sum_twist = 0.0;
    long num_positive_twist = 0;
    long i, j, nitem, nlen, ip = 0;
    FILE *fp;

    fp = open_file(inpfile, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", &pars->num_bp) != 1)
        fatal("error reading parameter file: number of base-pairs\n");
    if (pars->num_bp <= 0)
        fatal("number of base-pairs [%ld] <= 0\n", pars->num_bp);

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", &pars->is_helical) != 1)
        fatal("error reading parameter file: step or helical option\n");

    if (fgets(str, sizeof str, fp) == NULL)
        fatal("error reading parameter file: parameter title line\n");

    pars->bp_seq = cmatrix(1, pars->num_bp, 0, 2);
    pars->bp_par = dmatrix(1, pars->num_bp, 1, 6);
    pars->step_par = dmatrix(1, pars->num_bp, 1, 6);
    pars->pidx = lvector(1, pars->num_bp);

    for (i = 1; i <= pars->num_bp; i++) {
        if (fgets(str, sizeof str, fp) == NULL)
            fatal("error reading parameter file: parameters list <%ld: %s>\n", i, str);

        p = strrchr(str, '#');
        if (p)  /* remove in-line comment starting with '#' */
            *p = '\0';

        nitem = itemize(str, item, 13);
        nlen = strlen(item[0]);  /* keep case as is: g for +G etc */
        pars->bp_seq[i][1] = *item[0];
        if (nlen == 1)
            ++pars->is_single;
        else if (nlen == 3) {
            pars->bp_seq[i][0] = *(item[0] + 1);  /* - or + */
            if (pars->bp_seq[i][0] != '+')  /* default to anti-parallel */
                pars->bp_seq[i][0] = '-';
            pars->bp_seq[i][2] = *(item[0] + 2);
        } else
            fatal("error reading parameter file: base-pair format\n");

        if (nitem == 6 || nitem == 7) {  /* without bp parameters */
            if (pars->num_bp == 1) {  /* a single bp: taken as base-pair parameters */
                for (j = 1; j <= 6; j++)
                    if (sscanf(item[j], "%lf", &pars->bp_par[i][j]) != 1)
                        fatal("error reading base-pair parameters\n");
            } else {  /* taken as step parameters: bp_par = 0.0 */
                for (j = 1; j <= 6; j++)
                    if (sscanf(item[j], "%lf", &pars->step_par[i][j]) != 1)
                        fatal("error reading step/helical parameters\n");
            }
        } else if (nitem == 12 || nitem == 13) {  /* with bp parameters */
            for (j = 1; j <= 6; j++)
                if (sscanf(item[j], "%lf", &pars->bp_par[i][j]) != 1)
                    fatal("error reading base-pair parameters\n");
            for (j = 1; j <= 6; j++)
                if (sscanf(item[j + 6], "%lf", &pars->step_par[i][j]) != 1)
                    fatal("error reading step/helical parameters\n");
        } else
            fatal("wrong number of parameters per line\n");

        if (nitem == 7 || nitem == 13) {
            if (sscanf(item[nitem], "%ld", &pars->pidx[i]) != 1)
                fatal("error reading P xyz index\n");
        } else
            pars->pidx[i] = 2;  /* B-DNA (index 2) as the default */
    }

    close_file(fp);

    if (pars->is_single && pars->is_single != pars->num_bp)
        fatal("sorry, no mixed single helix and duplex allowed\n");

    /* set +x-axis direction */
    for (i = 1; i <= pars->num_bp; i++) {
        twist = pars->step_par[i][6];
        sum_twist += twist;
        if (twist >= 0)
            num_positive_twist++;
        if (pars->bp_seq[i][0] == '+')
            ip++;
    }
    pars->xdir = twist < 0.0 && !num_positive_twist;  /* as in 1xvr */
    if (ip == pars->num_bp)
        pars->parallel = 1;
}

int main(int argc, char *argv[])
{
    char BDIR[BUF512], inpfile[BUF512], outfile[BUF512];
    char user_input[BUF512];
    char *options[] = { "-ATOMIC", "-BASE_P", "-BLOCK1", "-BLOCK2" };
    long num_bp, num_atoms, num_max_per_residue, neg_xdir = 0, xml = 0;
    struct_pars pars;
    time_t time0;

    time(&time0);

    remove_file(REF_FILE);

    set_my_globals(argv[0]);

    strcpy(user_input, options[2]);  /* block1 */
    rebuild_cmdline(argc, argv, options, user_input, inpfile, outfile, &neg_xdir, &xml);

    init_pars(&pars);
    read_pars(inpfile, &pars);
    num_bp = pars.num_bp;

    if (pars.xdir && !neg_xdir)
        fprintf(stderr, "Twist < 0.0, taken as Z-DNA and reverse x- and z-axes\n");
    if (neg_xdir)
        pars.xdir = !pars.xdir;

    if (!strcmp(user_input, options[0])) {
        get_BDIR(BDIR, "Atomic_A.pdb");
        num_PDB_atoms(num_bp, pars.is_single, pars.bp_seq, BDIR, &num_atoms, &num_max_per_residue);
        if (pars.is_single)
            atomic_pdb1(num_bp, num_atoms, num_max_per_residue, pars.is_helical,
                        pars.xdir, pars.bp_seq, pars.step_par, BDIR, xml, outfile);
        else
            atomic_pdb2(pars.parallel, num_bp, num_atoms, num_max_per_residue,
                        pars.is_helical, pars.xdir, pars.bp_seq, pars.bp_par,
                        pars.step_par, BDIR, xml, outfile);
    } else if (!strcmp(user_input, options[1])) {
        get_BDIR(BDIR, "Atomic_A.pdb");
        if (pars.is_single)
            fatal("<base_p> option only valid for a duplex\n");
        else {
            num_PDB_atoms(num_bp, pars.is_single, pars.bp_seq, BDIR, &num_atoms,
                          &num_max_per_residue);
            atomic_base_p(pars.parallel, num_bp, num_atoms, num_max_per_residue,
                          pars.is_helical, pars.xdir, pars.bp_seq, pars.bp_par,
                          pars.step_par, BDIR, pars.pidx, xml, outfile);
        }
    } else if (!strcmp(user_input, options[2])) {
        get_BDIR(BDIR, "Block_R.alc");
        block_alc1(num_bp, pars.is_single, pars.is_helical, pars.xdir, pars.bp_seq,
                   pars.step_par, BDIR, outfile);
    } else {
        if (pars.is_single)
            fatal("only one block per base allowed in a single helix\n");
        else {
            get_BDIR(BDIR, "Block_R.alc");
            block_alc2(num_bp, pars.is_helical, pars.xdir, pars.bp_seq, pars.bp_par,
                       pars.step_par, BDIR, outfile);
        }
    }

    free_pars(&pars);

    clear_my_globals();

    print_used_time(time0);

    return 0;
}
