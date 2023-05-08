#include "x3dna.h"

typedef struct {
    char inp1[BUF512];
    char inp2[BUF512];
    char out[BUF512];
} struct_args;

static void comb_str_usage(void)
{
    help3dna_usage("comb_str");
}

static void comb_str_cmdline(int argc, char *argv[], struct_args * args)
{
    long i;

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        comb_str_usage();  /* no other options allowed */
    }

    if (argc == i + 3) {
        strcpy(args->inp1, argv[i]);
        strcpy(args->inp2, argv[i + 1]);
        strcpy(args->out, argv[i + 2]);
    } else
        comb_str_usage();
}

static void comb_alcfile(char *inpfile1, char *inpfile2, char *outfile)
{
    char **AtomName, **AtomName2;
    double **xyz, **xyz2;
    long i, j, k, nbond, nbond1, nbond2, num, num1, num2;
    long *ibase, *ibase2;
    long **linkage, **linkage2;

    /* get the number of atoms and bonds */
    get_alc_nums(inpfile1, &num1, &nbond1);
    get_alc_nums(inpfile2, &num2, &nbond2);
    num = num1 + num2;
    nbond = nbond1 + nbond2;

    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 1, 2);
    AtomName2 = cmatrix(1, num2, 0, 2);
    xyz2 = dmatrix(1, num2, 1, 3);
    ibase2 = lvector(1, num2);
    linkage2 = lmatrix(1, nbond2, 1, 2);

    read_alc(inpfile1, &num1, &nbond1, AtomName, xyz, ibase, linkage);
    read_alc(inpfile2, &num2, &nbond2, AtomName2, xyz2, ibase2, linkage2);
    for (i = 1; i <= num2; i++) {
        k = num1 + i;
        strcpy(AtomName[k], AtomName2[i]);
        ibase[k] = ibase2[i];
        cpxyz(xyz2[i], xyz[k]);
    }
    for (i = 1; i <= nbond2; i++)
        for (j = 1; j <= 2; j++)
            linkage[nbond1 + i][j] = num1 + linkage2[i][j];

    write_alc(num, nbond, AtomName, xyz, ibase, linkage, outfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_alc(num2, nbond2, AtomName2, xyz2, ibase2, 1, linkage2);
}

/* read in the CONECT records from a PDB file.
 * array `connect' has the following format:
 *     col# 1-6: maximum 6 connections
 *     col# 7: total number of connections
 * return value:
 *     0 if no such CONECT records or its # < the # of atom records
 *     1 otherwise return */
static long get_pdbconnect(long num, long **connect, char *pdbfile)
{
    char str[BUF512], *item[MBASES];
    long i, k, nitem, nlen, nc = 0;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        nlen = upperstr(str);
        if (nlen >= 6 && !strncmp(str, "CONECT", 6)) {  /* at least 6 characters */
            nitem = itemize(str + 6, item, MBASES);
            if (!nitem)
                fprintf(stderr, "Record <%s> has no contents\n", str);
            else if (nitem > 6)  /* maximum 6 linkages */
                nitem = 6;
            if (sscanf(item[0], "%ld", &k) != 1)
                fatal("error reading atom serial number in CONECT record\n");
            for (i = 1; i <= nitem; i++)
                if (sscanf(item[i], "%ld", &connect[k][i]) != 1)
                    fatal("error reading connected atom serial # in CONECT record\n");
            nc++;
            connect[k][7] = nitem;
        }
    }
    close_file(fp);

    return (!nc || nc < num) ? 0 : 1;
}

static void copy_pdb_coordinates(char *inppdb, long i0, FILE * fpo)
{
    char line[BUF512], snum[BUF32];
    long k = i0;
    FILE *fpi;

    fpi = open_file(inppdb, "r");
    while (fgets(line, sizeof line, fpi) != NULL) {
        if (str_pmatch(line, "ATOM  ") || str_pmatch(line, "HETATM")) {
            k++;
            sprintf(snum, "%5ld", k);
            strncpy(line + 6, snum, 5);
            fprintf(fpo, "%s", line);
        }
    }
    close_file(fpi);
}

static void comb_pdbfile(long num1, char *inpfile1, char *inpfile2, char *outfile)
{
    char *ChainID, **AtomName, **ResName, **Miscs;
    char *ChainID2, **AtomName2, **ResName2, **Miscs2;
    double **xyz, **xyz2;
    long i, j, k, num, num2;
    long *ResSeq, **connect;
    long *ResSeq2, **connect2;

    num2 = number_of_atoms(inpfile2, 1, "*");
    num = num1 + num2;
    connect = lmatrix(1, num, 1, 7);
    connect2 = lmatrix(1, num2, 1, 7);

    if (!get_pdbconnect(num1, connect, inpfile1) ||  /* using num1 */
        !get_pdbconnect(num2, connect2, inpfile2)) {
        FILE *fp;

        fp = open_file(outfile, "w");
        fprintf(fp, "NOTE  combined from '%s' and '%s'\n", inpfile1, inpfile2);
        copy_pdb_coordinates(inpfile1, 0, fp);
        fprintf(fp, "NOTE  -- '%s'\n", inpfile2);
        copy_pdb_coordinates(inpfile2, num1, fp);
        fprintf(fp, "END\n");

        close_file(fp);

    } else {
        AtomName = cmatrix(1, num, 0, 4);
        ResName = cmatrix(1, num, 0, 3);
        ChainID = cvector(1, num);
        ResSeq = lvector(1, num);
        xyz = dmatrix(1, num, 1, 3);
        Miscs = cmatrix(1, num, 0, NMISC);
        read_pdb(inpfile1, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");

        AtomName2 = cmatrix(1, num2, 0, 4);
        ResName2 = cmatrix(1, num2, 0, 3);
        ChainID2 = cvector(1, num2);
        ResSeq2 = lvector(1, num2);
        xyz2 = dmatrix(1, num2, 1, 3);
        Miscs2 = cmatrix(1, num2, 0, NMISC);
        read_pdb(inpfile2, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz2, Miscs2, 1, "*");

        for (i = 1; i <= num2; i++) {
            k = num1 + i;
            strcpy(AtomName[k], AtomName2[i]);
            strcpy(ResName[k], ResName2[i]);
            ChainID[k] = ChainID2[i];
            ResSeq[k] = ResSeq2[i];
            cpxyz(xyz2[i], xyz[k]);
            strcpy(Miscs[k], Miscs2[i]);
            connect[k][7] = connect2[i][7];
            for (j = 1; j <= connect2[i][7]; j++)
                connect[k][j] = num1 + connect2[i][j];
        }

        write_pdbcnt(num, AtomName, ResName, ChainID, ResSeq, xyz, connect, outfile);

        free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
        free_pdb(num2, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz2, Miscs2);
    }

    free_lmatrix(connect, 1, num, 1, 7);
    free_lmatrix(connect2, 1, num2, 1, 7);
}

/* combine two ALCHEMY or PDB files with CONECT record. If either PDB file has no
 *    CONECT record at all or not enough records, simply concatenate them as if a
 *    `cat' command is used */
int main(int argc, char *argv[])
{
    long num1;
    struct_args args;

    set_my_globals(argv[0]);

    comb_str_cmdline(argc, argv, &args);

    num1 = number_of_atoms(args.inp1, 1, "*");
    (num1) ? comb_pdbfile(num1, args.inp1, args.inp2, args.out) :
        comb_alcfile(args.inp1, args.inp2, args.out);

    clear_my_globals();

    return 0;
}
