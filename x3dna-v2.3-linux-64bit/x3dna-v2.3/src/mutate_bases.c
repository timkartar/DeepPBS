#include "x3dna.h"

char fixed_mutfile[BUF512] = "mutations.dat";

#define DUMCHAR '@'
#define DUMNAME "@@@"
typedef struct {
    char chainID;
    char resName[4];
    char iCode;
    long resSeq;

    char mutName[4];  /* 3-char */

    /* -1: required chainID or resSeq not specified
     * -2: duplicate
     * -3: residue name to be mutated not in "baselist.dat"
     * -4: residue for mutation is identical
     * -5: residue to be mutated not in PDB file */
    long serialNum;
    char mutBase;  /* 1-char */
} mutInfo;

typedef struct {
    char pdbfile[BUF512];
    char outfile[BUF512];
    char mutfile[BUF512];
    long list;
    long enum_base;
} struct_args;

static void mutate_bases_usage(void)
{
    help3dna_usage("mutate_bases");
}

static void set_defaults(struct_args * args)
{
    strcpy(args->pdbfile, "");
    strcpy(args->outfile, "");
    strcpy(args->mutfile, "");
    args->list = FALSE;
    args->enum_base = FALSE;
}

static void mutate_bases_cmdline(int argc, char *argv[], struct_args * args)
{
    char c;
    long i, j;

    if (argc < 4)
        mutate_bases_usage();

    set_defaults(args);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        if (str_pmatch(argv[i], "-e")) {
            args->enum_base = TRUE;
            continue;
        }

        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++) {  /* skip - */
            c = argv[i][j];
            if (c == 'L')
                args->list = TRUE;
            else if (c == 'E')
                args->enum_base = TRUE;
            else
                mutate_bases_usage();
        }
    }

    if (args->enum_base) {  /* just for enumeration */
        if (argc == i + 2) {
            strcpy(args->pdbfile, argv[i]);
            strcpy(args->outfile, argv[i + 1]);
        } else
            mutate_bases_usage();

        if (!exist_file(args->pdbfile))
            fatal("PDB file <%s> does not exist!\n", args->pdbfile);
        return;
    }

    if (argc == i + 3) {
        strcpy(args->mutfile, argv[i]);
        strcpy(args->pdbfile, argv[i + 1]);
        strcpy(args->outfile, argv[i + 2]);
    } else
        mutate_bases_usage();

    if (!exist_file(args->pdbfile))
        fatal("PDB file <%s> does not exist!\n", args->pdbfile);

    if (args->list && !exist_file(args->mutfile))
        fatal("Mutation file <%s> does not exist!\n", args->mutfile);
}

static void base_idmsg(char chain_id, long res_seq, char icode, char *rname, char *idmsg)
{
    char snum[BUF32], rname_cp[BUF32];
    long i;

    if (res_seq == LONG_MAX)
        strcpy(snum, "NONE");
    else {
        sprintf(snum, "%4ld", res_seq);
        for (i = 0; i < 4; i++)
            if (snum[i] == ' ')
                snum[i] = '.';
    }

    if (chain_id == ' ')
        chain_id = '-';

    if (icode == ' ')
        icode = '_';

    strcpy(rname_cp, rname);
    for (i = 0; i < 3; i++)
        if (rname_cp[i] == ' ')
            rname_cp[i] = '.';

    sprintf(idmsg, "%c:%4s%c:[%s]", chain_id, snum, icode, rname_cp);
}

static void print_mutInfo(long num_mut, mutInfo * mutation)
{
    char idmsg[BUF512];
    long i, num = 0;

    mutInfo *m;

    for (i = 1; i <= num_mut; i++) {
        m = &mutation[i];
        base_idmsg(m->chainID, m->resSeq, m->iCode, m->resName, idmsg);
        fprintf(stderr, "%4ld\t%s\t ===> [%s.%c]\t", i, idmsg, m->mutName, m->mutBase);
        if (m->serialNum < 0) {
            switch (m->serialNum) {
            case -1:
                fprintf(stderr, "required chainID/resnum/mutation not specified\n");
                break;
            case -2:
                fprintf(stderr, "duplicate of previous entry\n");
                break;
            case -3:
                fprintf(stderr, "residue name to be mutated not in 'baselist.dat'\n");
                break;
            case -4:
                fprintf(stderr, "identical residue for mutation\n");
                break;
            case -5:
                fprintf(stderr, "residue to be mutated not in the PDB file\n");
                break;
            default:
                fatal("how could it possible?\n");
            }
        } else {
            fprintf(stderr, "---done---\n");
            num++;
        }
    }
    fprintf(stderr, "    Number of mutations: %ld\n", num);
}

static void initialize_mutInfo_entry(mutInfo * m)
{
    m->chainID = DUMCHAR;
    strcpy(m->resName, DUMNAME);
    m->iCode = DUMCHAR;
    m->resSeq = LONG_MAX;
    strcpy(m->mutName, DUMNAME);

    m->serialNum = 0;
    m->mutBase = DUMCHAR;
}

static void initialize_mutInfo(long num, mutInfo * ids)
{
    long i;

    for (i = 1; i <= num; i++)
        initialize_mutInfo_entry(&ids[i]);
}

static mutInfo *allocate_memory_for_mutInfo(long num)
{
    mutInfo *ids;

    ids = (mutInfo *) malloc((num + 1) * sizeof(mutInfo));

    if (ids == NULL)
        fatal("malloc failure for allocate_memory_for_mutInfo(%ld)\n", num);

    initialize_mutInfo(num, ids);

    return ids;
}

static void free_mutInfo(mutInfo * ids)
{
    free(ids);
}

static void populate_mutation_entry(char *par_val, mutInfo * m)
{
    char *p, *items[BUF512], val[BUF512], str0[BUF512];
    long i, k, nitem, NITEM_MAX = 5;

    strcpy(str0, par_val);
    null_line_comment(par_val);
    upperstr(par_val);  /* case-insensitive */
    nitem = item_list(par_val, items, NITEM_MAX, WSPACES);

    if (nitem < 3) {
        m->serialNum = -1;
        fprintf(stderr, "%s\tmiss required fields: chain, seqnum or mutation\n", str0);
    }

    for (i = 1; i <= nitem; i++) {
        p = items[i];
        get_strvalue(p, val, FALSE);
        k = strlen(val);

        if (str_pmatch(p, "C")) {
            if (k == 0)
                fprintf(stderr, "no chain ID specified: <%s>\n", p);
            else if (k > 1)
                fprintf(stderr, "chain ID should be 1-char: <%s>\n", p);
            else {
                m->chainID = val[0];
                if (m->chainID == '_')
                    m->chainID = ' ';
            }

        } else if (str_pmatch(p, "N") || str_pmatch(p, "RESNA")) {
            if (k > 3)
                fprintf(stderr, "residue name more than 3-char: <%s>\n", p);
            else
                sprintf(m->resName, "%3s", val);

        } else if (str_pmatch(p, "I")) {
            if (k == 0)
                fprintf(stderr, "no insert code specified: <%s>\n", p);
            else if (k > 1)
                fprintf(stderr, "insert code should be 1-char: <%s>\n", p);
            else
                m->iCode = val[0];

        } else if (str_pmatch(p, "S") || str_pmatch(p, "RESS") || str_pmatch(p, "RESNU")) {
            m->resSeq = get_lvalue(p, LONG_MIN, LONG_MAX);

        } else if (str_pmatch(p, "M")) {
            if (k > 3)
                fprintf(stderr, "mutated residue name > 3-char: <%s>\n", p);
            else
                sprintf(m->mutName, "%3s", val);

        } else
            fprintf(stderr, "unrecognized option: <%s>\n", p);
    }
}

static void populate_mutation(char *mutfile, mutInfo * mutation)
{
    char *p0, *line;
    long num = 0;
    FILE *fp;

    fp = open_file(mutfile, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        line = ltrim(p0);

        if (is_skip_line(line)) {
            free(p0);
            continue;  /* skip empty and commented lines */
        }

        num++;
        populate_mutation_entry(line, &mutation[num]);

        free(p0);
    }

    close_file(fp);
}

static void check_mutation_required_fields(long num, mutInfo * mutation)
{
    long i;
    mutInfo *m;

    for (i = 1; i <= num; i++) {  /* check for required chainID and resSeq */
        m = &mutation[i];
        if (m->serialNum != 0)
            continue;
        if (m->chainID == DUMCHAR || m->resSeq == LONG_MAX || is_equal_string(m->mutName, DUMNAME)) {
            m->serialNum = -1;
            fprintf(stderr, "%ld\tmiss required fields: chain, seqnum or mutation\n", i);
        }
    }
}

static void check_mutation_duplicates(long num, mutInfo * mutation)
{
    char str[BUF512];
    long i, j;

    mutInfo *m, *n;

    for (i = 1; i < num; i++) {
        m = &mutation[i];
        if (m->serialNum != 0)
            continue;

        for (j = i + 1; j <= num; j++) {
            n = &mutation[j];
            if (n->serialNum != 0)
                continue;

            if ((m->chainID == n->chainID) && is_equal_string(m->resName, n->resName) &&
                (m->iCode == n->iCode) && (m->resSeq == n->resSeq)) {
                base_idmsg(m->chainID, m->resSeq, m->iCode, m->resName, str);
                fprintf(stderr, "%ld and %ld are duplicates: %s\n", i, j, str);
                m->serialNum = -2;
                fprintf(stderr, "%ld and %ld are duplicates\n", i, j);
            }
        }
    }
}

static void check_mutation_baselist(long num, mutInfo * mutation)
{
    char str[BUF512];
    long i, j;

    mutInfo *m;

    for (i = 1; i <= num; i++) {
        m = &mutation[i];
        if (m->serialNum != 0)
            continue;

        base_idmsg(m->chainID, m->resSeq, m->iCode, m->resName, str);
        for (j = 1; j <= Gvars.NUM_SBASE; j++)
            if (strncmp(Gvars.BASELIST[j], m->mutName, 3) == 0)
                break;

        if (j > Gvars.NUM_SBASE) {
            m->serialNum = -3;
            fprintf(stderr, "unrecognized mutation: %ld %s ---> *%s*\n", i, str, m->mutName);
        } else
            m->mutBase = Gvars.BASELIST[j][3];
    }
}

static void check_mutation_identical(long num, mutInfo * mutation)
{
    long i;

    mutInfo *m;

    for (i = 1; i <= num; i++) {
        m = &mutation[i];
        if (m->serialNum != 0)
            continue;
        if (is_equal_string(m->resName, m->mutName)) {
            m->serialNum = -4;
            fprintf(stderr, "%ld: mutation is the same as original\n", i);
        }
    }
}

static void check_mutation(long num, mutInfo * mutation)
{
    long i, num_ok = 0;

    check_mutation_required_fields(num, mutation);  /* -1 */
    check_mutation_duplicates(num, mutation);  /* -2 */
    check_mutation_baselist(num, mutation);  /* -3 */
    check_mutation_identical(num, mutation);  /* -4 */

    num_ok = 0;
    for (i = 1; i <= num; i++)
        if (mutation[i].serialNum == 0)
            num_ok++;

    if (num_ok != num) {
        fprintf(stderr, "%ld out of %ld mutation entries are valid.\n", num_ok, num);
        if (num_ok == 0)
            fatal("\tNO valid mutation entries -- check file '%s'!\n", fixed_mutfile);
    }
}

static long match_iCode(char mut_iCode, char res_iCode)
{
    if (mut_iCode == DUMCHAR)
        return TRUE;

    if (mut_iCode == res_iCode)
        return TRUE;

    return FALSE;
}

static long match_resName(char *mutName, char *resName)
{
    if (is_equal_case_string(mutName, DUMNAME))
        return TRUE;

    if (is_equal_case_string(mutName, resName))
        return TRUE;

    return FALSE;
}

static void match_mutInfo(long num_mut, mutInfo * mutation, long num_residue,
                          long **seidx, long *RY, char *ChainID, char **ResName,
                          char **Miscs, long *ResSeq)
{
    char iCode, str[BUF512], idmsg[BUF512];
    long i, ib, j, k;

    mutInfo *m;

    for (i = 1; i <= num_mut; i++) {
        m = &mutation[i];
        if (m->serialNum != 0)
            continue;
        base_idmsg(m->chainID, m->resSeq, m->iCode, m->resName, str);

        k = 0;
        for (j = 1; j <= num_residue; j++) {
            if (RY[j] < 0)
                continue;
            ib = seidx[j][1];
            iCode = Miscs[ib][2];
            if ((m->chainID == toupper((int) ChainID[ib])) && (m->resSeq == ResSeq[ib]) &&
                match_iCode(m->iCode, iCode) && match_resName(m->resName, ResName[ib])) {
                m->serialNum = j;
                k++;
            }
        }

        if (k == 0) {
            fprintf(stderr, "Mutation entry %ld %s has no PDB residue match\n", i, str);
            m->serialNum = -5;

        } else {
            ib = seidx[m->serialNum][1];
            if (k > 1) {
                fprintf(stderr, "Mutation entry %ld %s has %ld matches\n", i, str, k);
                base_idmsg(ChainID[ib], ResSeq[ib], Miscs[ib][2], ResName[ib], idmsg);
                fprintf(stderr, "\tlast match mutated: %s\n", idmsg);
            }
        }
    }
}

static void mutate_record(long ib, long ie, long *inum, char mutBase, char *mutName,
                          double *orien, double *org, char cid, long rnum,
                          char **AtomName, double **xyz, char **Miscs, FILE * fp)
{
    char rname[4], aname[5], str[BUF512];
    char spdb[BUF512], BDIR[BUF512], *sChainID, **sAtomName, **sResName;
    long i, j, is_dna, snum, *sResSeq;
    double txyz[4], **sxyz, **orien_cmtx;

    get_BDIR(BDIR, "Atomic_A.pdb");

    sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
    sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
    sChainID = cvector(1, NUM_RESIDUE_ATOMS);
    sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
    sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);

    if (!set_3letter_base_pdb(mutName, spdb))
        set_std_base_pdb(BDIR, FALSE, mutBase, spdb);
    snum = read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL, 1, "*");

    /* is_dna = Gvars.PDBV3 && is_dna_with_backbone(ib, ie, AtomName); */
    is_dna = FALSE;  /* respect user input */

    for (i = ib; i <= ie; i++) {  /* for backbone atoms, as is */
        if (is_baseatom(AtomName[i]))
            continue;

        deduce_misc(Miscs, AtomName, i, str);
        normalize_resName_atomName(is_dna, mutName, AtomName[i], rname, aname);

        fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n",
                (str[0] == 'A') ? "ATOM  " : "HETATM", ++*inum, aname,
                str[1], rname, cid, rnum, str[2], xyz[i][1], xyz[i][2], xyz[i][3], str + 3);
    }

    orien_cmtx = dmatrix(1, 3, 1, 3);
    orien2mst(orien, 0, orien_cmtx);

    for (i = 1; i <= snum; i++) {  /* for mutated base atom */
        if (!is_baseatom(sAtomName[i]))
            continue;

        for (j = 1; j <= 3; j++)
            txyz[j] = dot(sxyz[i], orien_cmtx[j]) + org[j];
        deduce_misc(NULL, sAtomName, i, str);
        normalize_resName_atomName(is_dna, mutName, sAtomName[i], rname, aname);

        if (is_equal_string(rname, "5CM") && is_equal_string(aname, " C7 "))
            strcpy(aname, " C5A");  /* PDB convention */

        fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n",
                (str[0] == 'A') ? "ATOM  " : "HETATM", ++*inum, aname, str[1],
                rname, cid, rnum, str[2], txyz[1], txyz[2], txyz[3], str + 3);
    }

    free_pdb(NUM_RESIDUE_ATOMS, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL);
    free_dmatrix(orien_cmtx, 1, 3, 1, 3);
}

static void mutate_pdb(long num_mut, mutInfo * mutation, long num_residue, long **seidx,
                       char **AtomName, char *ChainID, char **ResName, char **Miscs,
                       long *ResSeq, double **xyz, double **orien, double **org, char *outfile)
{
    char cid, idmsg[BUF512];
    long i, ib, ie, isMut, j, mok = 0, inum = 0, rnum;
    FILE *fp;

    mutInfo *m;

    fp = open_file(outfile, "w");
    fprintf(fp, "REMARK    %s\n", Gvars.X3DNA_VER);
    fprintf(fp, "REMARK    Mutated PDB File from 3DNA 'mutate_bases'\n");

    for (i = 1; i <= num_residue; i++) {
        isMut = FALSE;
        ib = seidx[i][1];
        ie = seidx[i][2];

        for (j = 1; j <= num_mut; j++) {
            m = &mutation[j];
            if (m->serialNum == i) {
                isMut = TRUE;
                break;
            }
        }

        if (isMut) {
            mok++;
            cid = ChainID[ib];
            rnum = ResSeq[ib];
            base_idmsg(m->chainID, m->resSeq, m->iCode, m->resName, idmsg);
            fprintf(fp, "REMARK    Mutation#%ld %s to [%s]\n", mok, idmsg, m->mutName);

            mutate_record(ib, ie, &inum, m->mutBase, m->mutName, orien[i], org[i], cid,
                          rnum, AtomName, xyz, Miscs, fp);

        } else
            pdb_record(seidx[i][1], seidx[i][2], &inum, 0, AtomName, ResName, ChainID,
                       ResSeq, xyz, Miscs, fp);
    }

    fprintf(fp, "END\n");

    close_file(fp);
}

static void list_bases(struct_args * args)
{
    char iCode, BDIR[BUF512], **nt_info;
    char *ChainID, *bseq, **AtomName, **ResName, **Miscs;
    double **xyz;
    long i, ib, num, hetatm = 1, num_residue;
    long *ResSeq, *RY, **seidx;
    FILE *fp;

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
    get_BDIR(BDIR, "Atomic_A.pdb");
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);

    nt_info = cmatrix(1, num_residue, 0, BUF32);
    populate_nt_info(num_residue, seidx, ResName, ChainID, ResSeq, Miscs, bseq, nt_info);

    fp = open_file(args->outfile, "w");
    fprintf(fp,
            "# add m=BASE_NAME (up-to three letters) to an entry for mutation\n"
            "#   e.g., change the line\n"
            "#           chain=A snum=2 name=C\n"
            "#         to\n"
            "#           chain=A snum=2 name=C mutation=G\n"
            "#   to mutate base C (on chain A and with residue number 2) to G\n\n");
    fprintf(fp, "# Empty or comment (starting with #s) lines are ignored\n\n");

    for (i = 1; i <= num_residue; i++) {
        if (RY[i] < 0)
            continue;
        ib = seidx[i][1];
        iCode = Miscs[ib][2];
        fprintf(fp, "chain=%c   snum=%-4ld  name=%-3s", ChainID[ib], ResSeq[ib],
                ltrim(ResName[ib]));
        if (iCode != ' ')
            fprintf(fp, " icode=%c", iCode);
        fprintf(fp, "  # %s %c\n", nt_info[i], bseq[i]);
    }
    close_file(fp);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cvector(bseq, 1, num_residue);
    free_lvector(RY, 1, num_residue);
    free_cmatrix(nt_info, 1, num_residue, 0, BUF32);
}

static void mutate_bases(struct_args * args)
{
    char *items[BUF512];
    char BDIR[BUF512], *ChainID, *bseq, **AtomName, **ResName, **Miscs;
    double **xyz, **orien, **org;
    long i, num, hetatm = 1, num_residue;
    long num_mut, nitem_read, nitem_max = 16;
    long *ResSeq, *RY, **seidx;
    FILE *fp;

    mutInfo *mutation;

    if (!args->list) {
        nitem_read = item_list(args->mutfile, items, nitem_max, "/;:");
        if (nitem_read == 0)
            fatal("%s contains no mutation entry\n", args->mutfile);

        fp = open_file(fixed_mutfile, "w");
        for (i = 1; i <= nitem_read; i++)
            fprintf(fp, "%s\n", items[i]);
        close_file(fp);

    } else
        strcpy(fixed_mutfile, args->mutfile);

    num_mut = get_line_number(fixed_mutfile, TRUE);
    mutation = allocate_memory_for_mutInfo(num_mut);
    populate_mutation(fixed_mutfile, mutation);

    check_mutation(num_mut, mutation);

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
    get_BDIR(BDIR, "Atomic_A.pdb");
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);

    match_mutInfo(num_mut, mutation, num_residue, seidx, RY, ChainID, ResName, Miscs, ResSeq);

    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 3);
    base_frame(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID, ResSeq, Miscs,
               xyz, BDIR, orien, org);

    print_mutInfo(num_mut, mutation);

    mutate_pdb(num_mut, mutation, num_residue, seidx, AtomName, ChainID, ResName,
               Miscs, ResSeq, xyz, orien, org, args->outfile);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cvector(bseq, 1, num_residue);
    free_lvector(RY, 1, num_residue);
    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 3);

    free_mutInfo(mutation);
}

int main(int argc, char *argv[])
{
    struct_args args;
    time_t time0;

    time(&time0);

    set_my_globals(argv[0]);

    mutate_bases_cmdline(argc, argv, &args);

    if (args.enum_base)
        list_bases(&args);
    else
        mutate_bases(&args);

    clear_my_globals();

    print_used_time(time0);

    return 0;
}
