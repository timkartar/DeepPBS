#include "x3dna.h"

typedef struct {
    long istart;
    long istep;
    long irna;  /* analyze using the RNA algorithm */
    long bz;  /* check for B-Z junction */
} struct_args;

/* process a structure */
static void process_cehs(char *inpfile, struct_args * args)
{
    char pdbfile[BUF512], outfile[BUF512];
    char *ChainID, **AtomName, **bp_seq, **ResName, **Miscs;
    double **xyz;
    long ip, ds, hetatm, num, num_bp, num_residue, bbexist = 0, parallel = 0;
    long *ResSeq;
    long **c6_c8, **chi, **pair_num, **phos, *RY;
    long **seidx, **sugar, **o3p_brk;
    FILE *fp;
    time_t run_time;

    /* get PDB file and pairing information */
    pair_num = read_input(inpfile, pdbfile, outfile, &ds, &num_bp, &ip, &hetatm);
    if (ds == 1 && !args->irna) {
        fprintf(stderr, "CEHS only applies to duplex structures\n");
        return;
    }

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

    /* base-pairing residue number checking */
    pair_checking(ip, ds, num_residue, pdbfile, &num_bp, pair_num);

    strcat(outfile, args->irna ? "r" : "c");  /* add c/r for CEHS/RNA parameters: *.outc/r */
    fp = open_file(outfile, "w");

    /* check strand direction: O3'-P linkage, parallel etc */
    o3p_brk = lmatrix(1, ds, 1, num_bp);
    drct_checking(ds, num_bp, pair_num, seidx, AtomName, xyz, &parallel, &bbexist, o3p_brk, fp);

    /* get base or base-pair sequence */
    bp_seq = cmatrix(0, ds, 1, num_bp);
    RY = lvector(1, num_residue);
    get_bpseq(ds, num_bp, pair_num, seidx, AtomName, ResName, ChainID, ResSeq, Miscs,
              xyz, bp_seq, RY);

    /* get atom list for each residue */
    phos = lmatrix(1, ds + 4, 1, num_bp);  /* also include O1P & O2P atoms */
    c6_c8 = lmatrix(1, ds, 1, num_bp);
    sugar = lmatrix(1, ds, 1, num_bp * 5);
    chi = lmatrix(1, ds, 1, num_bp * 4);
    atom_list(ds, num_bp, pair_num, seidx, RY, bp_seq, AtomName, ResName,
              ChainID, ResSeq, Miscs, phos, c6_c8, sugar, chi);

    fprintf(fp,
            args->irna ? "             Authentic RNA Parameters Based on Babcock & Olson"
            : "        Authentic CEHS/SCHNAaP Parameters [RC8-YC6 & Normal Vectors]");
    fprintf(fp, "\n    %s\n", Gvars.X3DNA_VER);
    print_sep(fp, '*', 76);
    fprintf(fp, "File name: %s\n", pdbfile);

    run_time = time(NULL);
    fprintf(fp, "Date and time: %s\n", ctime(&run_time));

    fprintf(fp, "Number of base-pairs: %ld\n", num_bp);
    fprintf(fp, "Number of atoms: %ld\n", num);

    print_sep(fp, '*', 76);
    print_pdb_title(pdbfile, "*", fp);

    if (args->irna) {  /* RNA */
        double **org, **orien;
        long *WC_info, str_type = 0;

        orien = dmatrix(1, ds, 1, num_bp * 9);
        org = dmatrix(1, ds, 1, num_bp * 3);
        WC_info = lvector(1, num_bp);

        ref_frames(ds, num_bp, pair_num, bp_seq, seidx, RY, AtomName, ResName,
                   ChainID, ResSeq, Miscs, xyz, fp, orien, org, WC_info, &str_type,
                   args->irna, o3p_brk);
        out_rna(ds, num_bp, bp_seq, pair_num[ds + 1], orien, org, fp);

        free_dmatrix(orien, 1, ds, 1, num_bp * 9);
        free_dmatrix(org, 1, ds, 1, num_bp * 3);
    } else {  /* CEHS/SCHNAaP */
        double *bp_orien, *bp_org;

        bp_orien = dvector(1, num_bp * 9);
        bp_org = dvector(1, num_bp * 3);
        cehs_pars(num_bp, args->istart, args->istep, pair_num, bp_seq, seidx, c6_c8, RY,
                  AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bp_orien, bp_org, args->bz, fp);

        /* SCHNAaP global helical parameters: unlikely to have > 1 helix */
        schnaap_global(num_bp, num, bp_seq, chi, xyz, bp_orien, bp_org, fp);

        free_dvector(bp_orien, 1, num_bp * 9);
        free_dvector(bp_org, 1, num_bp * 3);
    }

    close_file(fp);

    /* free allocated vectors & matrices */
    free_lmatrix(pair_num, 1, ds + 1, 1, num_bp);
    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_lmatrix(o3p_brk, 1, ds, 1, num_bp);
    free_cmatrix(bp_seq, 0, ds, 1, num_bp);
    free_lvector(RY, 1, num_residue);
    free_lmatrix(phos, 1, ds + 4, 1, num_bp);
    free_lmatrix(c6_c8, 1, ds, 1, num_bp);
    free_lmatrix(sugar, 1, ds, 1, num_bp * 5);
    free_lmatrix(chi, 1, ds, 1, num_bp * 4);
}

static void cehs_usage(void)
{
    help3dna_usage("cehs");
}

static void set_defaults(struct_args * args)
{
    args->istart = 1;
    args->istep = 1;
    args->irna = FALSE;
    args->bz = TRUE;  /* check for B-Z junction */
}

static void cehs_cmdline(int argc, char *argv[], struct_args * args)
{
    long i, j, k;

    set_defaults(args);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        if (lux_ncmatch(argv[i], "^--?bz")) {
            args->bz = set_switch_default_true(argv[i]);
            continue;
        }

        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'R')
                args->irna = TRUE;
            else if (argv[i][j] == 'S') {
                k = sscanf(argv[i], "-S=%ld,%ld", &args->istep, &args->istart);
                if (k == 2) {
                    if (args->istart < 0)
                        args->istart = -args->istart;  /* make it positive */
                } else if (k == 1)
                    args->istart = 1;  /* redundant: already initialized */
                else {
                    fprintf(stderr, "wrong format for setting step\n");
                    cehs_usage();
                }
                fprintf(stderr, "***start at %ld, with step size: %ld***\n", args->istart,
                        args->istep);
                break;  /* -s=istep,istart not combined with others */
            } else
                cehs_usage();
    }

    if (argc == i)
        cehs_usage();

    for (j = i; j < argc; j++) {
        if (strcmp(argv[j], "tmpfile"))  /* analyze multiple structures */
            fprintf(stderr, "\n......Processing structure #%ld: <%s>......\n", j - i + 1, argv[j]);
        process_cehs(argv[j], args);
    }
}

int main(int argc, char *argv[])
{
    struct_args args;
    time_t time0;

    time(&time0);

    set_my_globals(argv[0]);
    cehs_cmdline(argc, argv, &args);
    clear_my_globals();

    print_used_time(time0);

    return 0;
}
