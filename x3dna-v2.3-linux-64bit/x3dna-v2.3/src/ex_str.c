#include "x3dna.h"

typedef struct {
    char inpfile[BUF512];
    char outfile[BUF512];
    long inum;
    long NMR;
    long BIOU;
    long DELH;
} struct_args;

static void exstr_usage(void)
{
    help3dna_usage("ex_str");
}

#if 0
/* check if "pdbfile" is an NMR file based on EXPDTA record */
static void check_ifNMR(char *pdbfile)
{
    char str[BUF512];
    long nmr = 0;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        upperstr(str);
        if (!strncmp(str, "ATOM", 4) || !strncmp(str, "HETATM", 6) || !strncmp(str, "END", 3))
            break;
        if (!strncmp(str, "EXPDTA", 6) && strstr(str, "NMR")) {
            nmr = 1;
            break;
        }
    }
    if (!nmr)
        fprintf(stderr, "Warning: ***unlikely an authentic NMR file***\n");
}
#endif

static long get_1st_model(char *pdbfile)
{
    char str[BUF512];
    long inum = 1;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        upperstr(str);
        if (!strncmp(str, "MODEL ", 6) && (sscanf(str, "%*s %ld", &inum) == 1))
            break;
    }
    close_file(fp);

    return inum;
}

/* best representative NMR conformer, use #1 if not assigned in REMARK 210.
 * typical cases:
 *    BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE : NULL [pdb28sp]
 *    BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE : 1    [pdb2bj2]
 *    BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE : 3    [pdb2lef]
 *    BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE : 5    [pdb3php]
 *
 * for structures pdb1qby, pdb1djd & pdb1bjd, the REMARK 210 says the
 *    best representative conformer is #12 or #13, but the data file
 *    actually does not has that model. Reset inum to 1 for such cases. */
static long best_conformer(char *pdbfile)
{
    char str[BUF512], numstr[BUF512];
    long bestrec = 0, inum;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        upperstr(str);
        if (!strncmp(str, "ATOM", 4) || !strncmp(str, "HETATM", 6) || !strncmp(str, "END", 3))
            break;
        if (!strncmp(str, "REMARK 210 BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE :", 59)) {
            bestrec = 1;  /* has BEST REP ... record */
            if (sscanf(str + 59, "%ld", &inum) != 1) {
                fprintf(stderr, "Best representative conformer NULL: set to 1\n");
                inum = 1;
            } else {  /* Check if model inum really exists */
                sprintf(numstr, "%ld", inum);
                while (fgets(str, sizeof str, fp) != NULL)
                    if (!strncmp(str, "MODEL ", 6) && strstr(str, numstr) != NULL) {
                        bestrec = 2;
                        break;
                    }
                if (bestrec != 2) {
                    fprintf(stderr, "Best representative conformer #%ld does NOT exist"
                            ": set to 1\n", inum);
                    inum = 1;
                }
            }
            break;
        }
    }
    close_file(fp);

    if (!bestrec) {
        fprintf(stderr, "No best representative conformer assigned: set to #1\n");
        inum = get_1st_model(pdbfile);
    }

    return inum;
}

static void set_defaults(struct_args * args)
{
    strcpy(args->inpfile, "");
    strcpy(args->outfile, "");
    args->inum = 1;
    args->NMR = 0;
    args->BIOU = FALSE;
    args->DELH = FALSE;
}

static void exstr_cmdline(int argc, char *argv[], struct_args * args)
{
    long i;

    if (argc < 3)
        exstr_usage();

    set_defaults(args);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        upperstr(argv[i]);
        if (sscanf(argv[i], "-%ld", &args->inum) == 1) {  /* -6 */
            if (args->inum <= 0)  /* just output information, do not reset it */
                fprintf(stderr, "\t[i] model #: %ld [<= 0]\n", args->inum);
        } else if (!strcmp(argv[i], "-NMR"))
            args->NMR = 1;
        else if (!strncmp(argv[i], "-NMRB", 5))
            args->NMR = 2;
        else if (!strncmp(argv[i], "-BIOU", 5))
            args->BIOU = TRUE;
        else if (!strncmp(argv[i], "-DELH", 5))
            args->DELH = 1;
        else
            exstr_usage();
    }
    if (argc == i + 2) {
        strcpy(args->inpfile, argv[i]);
        strcpy(args->outfile, argv[i + 1]);
    } else
        exstr_usage();
#if 0
    if (args->NMR)
        check_ifNMR(args->inpfile);
#endif
    if (args->NMR == 2)
        args->inum = best_conformer(args->inpfile);
}

static void extract_model(char *inpfile, char *outfile, long inum)
{
    char str[BUF512];
    long found_str = 0, num0;
    FILE *fpi, *fpo;

    fpi = open_file(inpfile, "r");
    fpo = open_file(outfile, "w");

    while (fgets(str, sizeof str, fpi) != NULL) {
        if (!strncmp(str, "MODEL ", 6) && (sscanf(str, "%*s %ld", &num0) == 1) && num0 == inum) {
            print_pdb_title(inpfile, "*", fpo);
            found_str = 1;
            fprintf(stderr, "Model selected #%ld\n", inum);
            if (fputs(str, fpo) == EOF)
                fatal("error in writing selected structure\n");
            while (1) {
                if (fgets(str, sizeof str, fpi) != NULL) {
                    fputs(str, fpo);
                    if (!strncmp(str, "ENDMDL", 6)) {
                        found_str = 2;
                        goto FINISHED;
                    }
                } else if (feof(fpi))
                    goto FINISHED;
                else
                    fatal("error reading PDB data file\n");
            }
        }
    }

  FINISHED:
    if (!found_str) {
        fprintf(stderr, "Model #: %ld not found in file %s\n", inum, inpfile);
        if (inum == 1) {  /* only for the 1st unit */
            fprintf(stderr, "    Copy the whole file: %s ===> %s\n", inpfile, outfile);
            cpcat_file(inpfile, outfile, "copy");
        }
    } else if (found_str != 2) {
        fprintf(stderr, "Model #: %ld has no ENDMDL record in file %s\n", inum, inpfile);
        fputs("ENDMDL\n", fpo);
        fputs("END\n", fpo);
    }

    close_file(fpi);
    close_file(fpo);
}

#if 0
/* extract a structure unit from a 3DNA output */
static void extract_3DNA(char *inpfile, char *outfile, long inum)
{
    char str[BUF512], numstr[BUF512];
    long found_str = 0;
    FILE *fpi, *fpo;

    fpi = open_file(inpfile, "r");
    fpo = open_file(outfile, "w");

    sprintf(numstr, "Section #%4.4ld", inum);
    while (fgets(str, sizeof str, fpi) != NULL) {
        if (strstr(str, numstr) != NULL) {
            found_str = 1;
            fprintf(stderr, "Structure selected #%ld\n", inum);
            if (fputs(str, fpo) == EOF)
                fatal("error in writing selected structure\n");
            while (1) {
                if (fgets(str, sizeof str, fpi) != NULL) {
                    fputs(str, fpo);
                    if (!strncmp(str, "END", 3)) {
                        found_str = 2;
                        goto FINISHED;
                    }
                } else if (feof(fpi))
                    goto FINISHED;
                else
                    fatal("error reading PDB data file\n");
            }
        }
    }

  FINISHED:
    if (!found_str) {
        fprintf(stderr, "3DNA structure #: %ld not found in file %s\n", inum, inpfile);
        if (inum == 1) {  /* only for the 1st unit */
            fprintf(stderr, "    Copy the whole file: %s ===> %s\n", inpfile, outfile);
            cpcat_file(inpfile, outfile, "copy");
        }
    } else if (found_str != 2) {
        fprintf(stderr, "3DNA structure #: %ld has no END record in file %s\n", inum, inpfile);
        fputs("END\n", fpo);
    }

    close_file(fpi);
    close_file(fpo);
}
#endif

/* chain id for each biounit */
static long biounit_cid(char *pdbfile, char **chlist)
{
    char str[BUF512], *item[MAXCH];
    long i, n, num_bu = 0;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        upperstr(str);
        if (!strncmp(str, "ATOM", 4) || !strncmp(str, "HETATM", 6) || !strncmp(str, "END", 3))
            break;
        if (!strncmp(str, "REMARK 350 APPLY THE FOLLOWING TO CHAINS:", 41)) {
            if (++num_bu > BUF512)
                fatal("Number of biological units exceed upper limit\n");
            n = itemize(str + 41, item, MAXCH) + 1;
            for (i = 0; i < n; i++)
                chlist[num_bu][i] = item[i][0];
            chlist[num_bu][i] = '\0';
        }
    }
    close_file(fp);

    return num_bu;
}

/* extract a biological unit from an x-ray crystal structure, based on
 * REMARK 350 record */
static void extract_BIOU(char *inpfile, char *outfile, long inum)
{
    char str[BUF512], *ChainID;
    char **AtomName, **ResName, **Miscs, **chlist;
    double **xyz;
    long j, k, nb, num, *ResSeq;
    FILE *fp;

    /* read in the PDB file */
    num = number_of_atoms(inpfile, 1, "*");
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");

    chlist = cmatrix(1, BUF512, 0, MAXCH);
    nb = biounit_cid(inpfile, chlist);

    if (!nb) {
        fprintf(stderr, "File %s does NOT contain `REMARK 350' records\n", inpfile);
        fprintf(stderr, "    Copy the whole file: %s ===> %s\n", inpfile, outfile);
        cpcat_file(inpfile, outfile, "copy");
    } else {
        fprintf(stderr, "Number of biological units: %ld\n", nb);
        if (inum > nb) {
            fprintf(stderr, "%ld > total # of biological unit: %ld\n", inum, nb);
            fprintf(stderr, "    Please re-select your biological unit #\n");
        } else {
            fprintf(stderr, "Biological unit selected #%ld\n", inum);
            fp = open_file(outfile, "w");
            print_pdb_title(inpfile, chlist[inum], fp);
            fprintf(fp, "REMARK 350 BIOMOLECULE: %ld\n", inum);
            k = 0;
            for (j = 1; j <= num; j++)
                if (strchr(chlist[inum], ChainID[j])) {
                    strcpy(str, Miscs[j]);
                    fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n",
                            (str[0] == 'A') ? "ATOM  " : "HETATM", ++k, AtomName[j],
                            str[1], ResName[j], ChainID[j], ResSeq[j], str[2], xyz[j][1],
                            xyz[j][2], xyz[j][3], str + 3);
                }
            fprintf(fp, "END\n");
            close_file(fp);
        }
    }
    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_cmatrix(chlist, 1, BUF512, 0, MAXCH);
}

/* extract a structure from a multiple-structure PDB file (NMR, stacking.pdb etc) */
int main(int argc, char *argv[])
{
    struct_args args;

    set_my_globals(argv[0]);

    exstr_cmdline(argc, argv, &args);
    remove_file(args.outfile);

    if (args.BIOU)
        extract_BIOU(args.inpfile, args.outfile, args.inum);
    else
        extract_model(args.inpfile, args.outfile, args.inum);

    if (args.DELH) {
        cpcat_file(args.outfile, TMP_FILE, "copy");
        delH_pdbfile(TMP_FILE, args.outfile);
    }

    clear_my_globals();

    return 0;
}
