#include "x3dna.h"

/* global variables definition */
struct_Gvars Gvars;

static char *CMARKERS = "||x+0";  /* helix begin, middle, end, isolated bp, last */
static char *REBUILD_CIDS = "ABCDEFGHI";

static char *ptable[] = {
    " H", "HE",
    "LI", "BE", " B", " C", " N", " O", " F", "NE",
    "NA", "MG", "AL", "SI", " P", " S", "CL", "AR",
    " K", "CA", "SC", "TI", " V", "CR", "MN", "FE", "CO",
    "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
    "RB", "SR", " Y", "ZR", "NB", "MO", "TC", "RU", "RH",
    "PD", "AG", "CD", "IN", "SN", "SB", "TE", " I", "XE",
    "CS", "BA", "LA", "HF", "TA", " W", "RE", "OS", "IR",
    "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN",
    "FR", "RA", "AC", "RF", "DB", "SG", "BH", "HS", "MT",
    "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB",
    "DY", "HO", "ER", "TM", "YB", "LU",
    "TH", "PA", " U", "NP", "PU", "AM", "CM", "BK",
    "CF", "ES", "FM", "MD", "NO", "LR"
};

static void get_3dna_homedir(char *homedir)
{
    char *temp;

    if ((temp = getenv("X3DNA")) != NULL) {
        strcpy(homedir, temp);
        check_slash(homedir);  /* with ending slash */
        return;
    }

    if ((temp = getenv("HOME")) != NULL ||  /* Unix */
        (temp = getenv("HOMEDRIVE")) != NULL) {  /* PC */
        char fname[BUF512];
        FILE *fp;

        strcpy(homedir, temp);
        check_slash(homedir);
        strcat(homedir, "x3dna-v2.3/");
        sprintf(fname, "%sconfig/help3dna.dat", homedir);
        fp = fopen(fname, "r");
        if (fp != NULL) {
            close_file(fp);
            return;
        }
    }

    fprintf(stderr, "Please set the X3DNA environment variable\n");
    fatal("Run 'x3dna_setup', and visit http://forum.x3dna.org/ for more info.\n");
}

static void get_3dna_version(char *homedir, char *version)
{
    char *p0, *line, str[BUF512];
    FILE *fp;

    sprintf(str, "%sconfig/version", homedir);
    fp = open_file(str, "r");
    while ((p0 = my_getline(fp)) != NULL) {
        line = trim(p0);  /* keep the original address of p0 */
        if (!is_skip_line(line)) {  /* just the first line */
            strcpy(version, line);
            free(p0);
            break;
        }
        free(p0);
    }
    close_file(fp);
}

long string_contains_only_those_characters(char *str, char *chars_set)
{
    return (strspn(str, chars_set) == strlen(str)) ? TRUE : FALSE;
}

/* bpid: +2 means a Watson-Crick base-pair;
 * bpid > 0 (+1): means WC geometry, as in G-U wobble
 * zdir < 0: '-' otherwise, '+' */
void bpid_wc_str(long bpid, double zdir, char *wc)
{
    sprintf(wc, "%c%c%c", (bpid == 2) ? '-' : '*', (bpid > 0) ? '-' : '*',
            (zdir < 0.0) ? '-' : '+');
}

int case_strcmp(const char *s1, const char *s2)
{
    int i, c1, c2;

    for (i = 0; c1 = toupper((int) s1[i]), c2 = toupper((int) s2[i]), c1 == c2; i++)
        if (c1 == '\0')
            return 0;
    return c1 - c2;
}

int case_strncmp(const char *s1, const char *s2, long n)
{
    int i, c1, c2;

    for (i = 0; (c1 = toupper((int) s1[i]), c2 = toupper((int) s2[i]), c1 == c2) && i < n; i++)
        if (c1 == '\0')
            return 0;

    return (i >= n) ? 0 : c1 - c2;
}

char *case_strstr(const char *haystack, const char *needle)
{
    char *haystack_cp, *needle_cp, *p;

    haystack_cp = my_strdup(haystack);
    needle_cp = my_strdup(needle);

    lowerstr(haystack_cp);
    lowerstr(needle_cp);

    p = strstr(haystack_cp, needle_cp);

    if (p != NULL)
        p = (p - haystack_cp) + (char *) haystack;

    free_cvector(haystack_cp, 0, -1);
    free_cvector(needle_cp, 0, -1);

    return p;
}

char *case_strchr(const char *s, int c)
{
    char *str, *p;

    str = my_strdup(s);
    lowerstr(str);

    p = strchr(str, tolower(c));

    if (p != NULL)
        p = (p - str) + (char *) s;

    free_cvector(str, 0, -1);

    return p;
}

void kbd_input(char *msg)
{
    if (msg != NULL)
        fprintf(stderr, "%s", msg);
    fprintf(stderr, "Press return to continue ...\n");
    getchar();
}

long get_line_number(char *filename, long skips)
{
    char *p0, *line;
    long num = 0;
    FILE *fp;

    fp = open_file(filename, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        if (skips) {
            line = trim(p0);  /* keep the original address of p0 */
            if (is_skip_line(line)) {
                free(p0);
                continue;
            }
        }
        num++;
        free(p0);
    }

    close_file(fp);

    if (!num)
        fatal("File <%s> contains 0 valid records\n", filename);

    return num;
}

long is_empty_string(const char *str)
{
    return (str == NULL || strcmp(str, "") == 0);
}

long is_equal_string(const char *str1, const char *str2)
{
    return (strcmp(str1, str2) == 0);
}

long is_equal_case_string(const char *str1, const char *str2)
{
    return (case_strcmp(str1, str2) == 0);
}

void null_line_comment(char *str)
{
    char *p;

    p = strrchr(str, '#');
    if (p) {
        if (!isspace((int) *(p - 1)))
            fatal("wrong format [%s]: requiring a white space before #\n", str);
        *p = '\0';
    }
}

long is_comment_line(char *line)
{
    char *p = ltrim(line);

    if (*p == '#')
        return TRUE;
    else
        return FALSE;
}

long is_empty_line(char *line)
{
    char *p = trim(line);

    if (*p == '\0')
        return TRUE;
    else
        return FALSE;
}

/* no processing of line[] here */
long is_skip_line(char *line)
{
    if (strchr(SKIPS, *line))  /* line starting with # or empty */
        return TRUE;
    else
        return FALSE;
}

/* get a new name with base name from 'src' and extension 'ext' */
void bname_ext(char *src, char *ext, char *dst)
{
    char str[BUF512];

    bname_noext(src, str);
    sprintf(dst, "%s.%s", str, ext);
}

double get_point2line_perp_distance(double *pnt, double *line_p1, double *line_p2)
{
    long i;
    double d, p1_p2[4], p1_pnt[4];

    ddxyz(line_p1, line_p2, p1_p2);
    vec_norm(p1_p2);

    ddxyz(line_p1, pnt, p1_pnt);
    d = dot(p1_pnt, p1_p2);

    for (i = 1; i <= 3; i++)
        p1_pnt[i] -= d * p1_p2[i];

    return veclen(p1_pnt);
}

void get_tag_string_pair(char *prefix, char *tag, char *btag, char *etag)
{
    sprintf(btag, "<%s%s", prefix, tag);
    sprintf(etag, "</%s%s", prefix, tag);
}

void get_xml_tag(FILE * fpxml, char *prefix, char *line, char *connector, char *tag, char *tag_str)
{
    char *p0, *p1, *ep, btag[BUF512], etag[BUF512], temp[BUFBIG];
    long isok = FALSE;

    get_tag_string_pair(prefix, tag, btag, etag);

    if (!tag_match(prefix, line, tag))
        return;

    strcpy(tag_str, "");  /* initialize it */
    ep = strstr(line, etag);

    if (ep) {  /* tag contained in one line */
        isok = TRUE;
        *ep = '\0';
        strcpy(temp, strchr(line, '>') + 1);
        strcpy(tag_str, trim(temp));

    } else {  /* in multiple lines */
        strcpy(tag_str, strchr(line, '>') + 1);
        while ((p0 = my_getline(fpxml)) != NULL) {
            p1 = trim(p0);
            ep = strstr(p1, etag);
            if (ep) {  /* tag contained in one line */
                isok = TRUE;
                *ep = '\0';
                strcat(tag_str, connector);
                strcat(tag_str, p1);
                free(p0);
                break;
            } else {
                strcat(tag_str, connector);
                strcat(tag_str, p1);
            }
            free(p0);
        }
    }

    if (!isok)
        fatal("tag <%s> has end match\n", tag);
}

void print_xml_tag(FILE * fpxml, char *prefix, char *line, char *tag, char *otag, FILE * fp)
{
    char tag_str[BUFBIG];

    if (tag_match(prefix, line, tag)) {
        get_xml_tag(fpxml, prefix, line, " ", tag, tag_str);
        fprintf(fp, "<%s>%s</%s>\n", otag, tag_str, otag);
    }
}

void get_xml_tag_long(FILE * fpxml, char *prefix, char *line, char *tag, long *lval)
{
    char tag_str[BUFBIG];

    if (tag_match(prefix, line, tag)) {
        get_xml_tag(fpxml, prefix, line, " ", tag, tag_str);
        *lval = cvt2long(tag_str);
    }
}

void get_xml_tag_double(FILE * fpxml, char *prefix, char *line, char *tag, double *dval)
{
    char tag_str[BUFBIG];

    if (tag_match(prefix, line, tag)) {
        get_xml_tag(fpxml, prefix, line, " ", tag, tag_str);
        *dval = cvt2double(tag_str);
    }
}

long tag_match(char *prefix, char *line, char *tag)
{
    char btag[BUF512];

    sprintf(btag, "<%s%s", prefix, tag);

    return str_pmatch(line, btag) && !strstr(line, "/>");  /* not short-handed form */
}

void extract_attribute(char *line, char *attr, char *attr_val, long to_lower)
{
    char *p0, *p1, *p2, str[BUF512];

    sprintf(str, "%s=\"", attr);

    p0 = strstr(line, str);
    if (!p0)
        fatal("line <%s> does not contain attribute: \"%s\"\n", line, attr);

    p1 = p0 + strlen(str);

    p2 = strchr(p1, '"');
    if (!p2)
        fatal("line <%s> has no matching \"\n", line);
    *p2 = '\0';

    strcpy(attr_val, p1);

    if (to_lower)
        lowerstr(attr_val);
}

void print_ptable()
{
    char str[BUF512];
    long i;

    for (i = 0; i <= Gvars.NUM_ELE; i++) {
        sprintf(str, "%2.2s..", Gvars.ATOM_NAMES[i]);
        if (str[0] == ' ')
            str[0] = '.';
        fprintf(stderr, "%4.4s  %2.2s   # %3ld\n", str, Gvars.ATOM_NAMES[i], i + 1);
    }
}

long set_switch_default_true(char *option)
{
    char str[BUF512];
    long switch_set = TRUE;

    if (!strchr(option, '='))  /* just a switch */
        return switch_set;

    get_strvalue(option, str, FALSE);
    lowerstr(str);

    if (str_pmatch(str, "of") || str_pmatch(str, "0") ||  /* off | 0 | no | false */
        str_pmatch(str, "n") || str_pmatch(str, "f"))
        switch_set = FALSE;

    return switch_set;
}

static void check_required_option(char *option, char *invalid_str, char *msg)
{
    if (is_equal_string(option, invalid_str)) {
        fprintf(stderr, "missing required option: %s\n", msg);
        fatal("Please try \"%s -h\" for usage information\n", Gvars.PROGNAME);
    }
}

void check_required_file(char *filename, char *invalid_str, char *msg)
{
    check_required_option(filename, invalid_str, msg);

    if (!exist_file(filename))
        fatal("required file <%s> does NOT exist!\n", filename);
}

void define_frame_by_3atoms(long num, double **xyz, double **refmat)
{
    double d, tmpx[4], tmpy[4], tmpz[4];
    long i;

    if (num < 3)
        return;

    cpxyz(xyz[1], refmat[4]);  /* 1st atom as origin */
    ddxyz(xyz[1], xyz[2], tmpx);  /* 1 ---> 2 as x-axis */

    for (i = 3; i <= num; i++) {
        ddxyz(xyz[1], xyz[i], tmpy);
        d = magang(tmpx, tmpy);
        if (dval_in_range(d, 5, 175))  /* nonlinear */
            break;
    }

    if (i > num) {
        fprintf(stderr, "too few non-linear atoms to define orientation\n");
        return;
    }

    vec_norm(tmpx);
    vec_orth(tmpy, tmpx);  /* orthogonal component as y-axis */
    cross(tmpx, tmpy, tmpz);

    x_y_z_2_mtx(tmpx, tmpy, tmpz, refmat);  /* column vector for x, y & z */
}

/* get a match of partial argument: e.g., -xml=2 */
long get_xmlArgNumber(char *str, char *pstr)
{
    long lval = 0;

    if (str_pmatch(str, pstr)) {
        lval = 1;
        if (strchr(str, '=') != NULL)
            lval = 2;
    }

    return lval;
}

/* get the local reference frame for each base. only the ring atoms are
   included in least-squares fitting. similar to ref_frames() */
void base_frame(long num_residue, char *bseq, long **seidx, long *res_type,
                char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                char **Miscs, double **xyz, char *BDIR, double **orien, double **org)
{
    static char *RingAtom[] = { RA_LIST };
    char idmsg[BUF512], sidmsg[BUF512], spdb[BUF512];
    char *sChainID, **sAtomName, **sResName, **sMiscs;
    double **sxyz, **eRing_xyz, **sRing_xyz, **fitted_xyz, **R;
    long i, ib, ie, j, RingAtom_num, RR9 = 9;
    long exp_katom, std_katom, nmatch, snum, *sResSeq;

    eRing_xyz = dmatrix(1, RR9, 1, 3);
    sRing_xyz = dmatrix(1, RR9, 1, 3);

    sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
    sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
    sChainID = cvector(1, NUM_RESIDUE_ATOMS);
    sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
    sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
    sMiscs = cmatrix(1, NUM_RESIDUE_ATOMS, 0, NMISC);

    fitted_xyz = dmatrix(1, RR9, 1, 3);
    R = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= num_residue; i++) {
        if (res_type[i] < 0)
            continue;  /* non-bases */

        ib = seidx[i][1];
        ie = seidx[i][2];
        get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);

        RingAtom_num = (res_type[i] == 1) ? 9 : 6;
        set_std_base_pdb(BDIR, FALSE, bseq[i], spdb);
        snum = read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, sMiscs, 1, "*");
        sprintf(sidmsg, "in standard base: %s", spdb);

        nmatch = 0;
        for (j = 0; j < RingAtom_num; j++) {
            exp_katom = find_1st_atom(RingAtom[j], AtomName, ib, ie, idmsg);
            std_katom = find_1st_atom(RingAtom[j], sAtomName, 1, snum, sidmsg);
            if (exp_katom && std_katom) {
                ++nmatch;
                cpxyz(xyz[exp_katom], eRing_xyz[nmatch]);
                cpxyz(sxyz[std_katom], sRing_xyz[nmatch]);
            }
        }

        (void) ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, org[i]);
        mst2orien(orien[i], 0, R);
    }

    free_pdb(NUM_RESIDUE_ATOMS, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, sMiscs);
    free_dmatrix(eRing_xyz, 1, RR9, 1, 3);
    free_dmatrix(sRing_xyz, 1, RR9, 1, 3);
    free_dmatrix(fitted_xyz, 1, RR9, 1, 3);
    free_dmatrix(R, 1, 3, 1, 3);
}

/* get residue index, color index and main chain atoms for the peptide bond */
void peptide_info(long num_residue, long **seidx, char **AtomName, char **ResName,
                  char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                  long *res_type, long *cidx, long *mchain)
{
    static char *SAA[] = { AA_LIST };
    char *col_code = "GGGGGGGAAAATTTTTTTTCX";
    char *cmn_base = CX_LIST, idmsg[BUF512];
    long i, ib, ie, ioffset, j, num_saa;

    num_saa = sizeof SAA / sizeof SAA[0] - 1;
    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];
        res_type[i] = residue_ident(AtomName, xyz, Miscs, ib, ie);
        if (res_type[i] != -1)  /* non-amino-acid */
            continue;
        get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);
        ioffset = (i - 1) * 4;
        mchain[ioffset + 1] = find_1st_atom(" N  ", AtomName, ib, ie, idmsg);
        mchain[ioffset + 2] = find_1st_atom(" CA ", AtomName, ib, ie, idmsg);
        mchain[ioffset + 3] = find_1st_atom(" C  ", AtomName, ib, ie, idmsg);
        mchain[ioffset + 4] = find_1st_atom(" O  ", AtomName, ib, ie, idmsg);
        for (j = 0; j <= num_saa; j++)
            if (!strcmp(ResName[ib], SAA[j]))
                break;
        cidx[i] = strchr(cmn_base, col_code[j]) - cmn_base;
    }
}

static void align_block_xyz(double *orien_i, double *org_i, double *xyz_j, double *new_xyz)
{
    long i, j;
    double temp[4];

    for (i = 1; i <= 3; i++) {
        for (j = 1; j <= 3; j++)
            temp[j] = orien_i[i + (j - 1) * 3];
        new_xyz[i] = dot(xyz_j, temp) + org_i[i];
    }
}

void base_blks(long num_residue, long *res_type, double **orien, double **org,
               char *bseq, char *BDIR, char *alcfile)
{
    char filename[BUF512], **AtomName, **tAtomName;
    double **txyz, **xyz;
    long bidx, i, j, k, nbond, num, tnbond, tnum;
    long ia = 0, ib = 0, num_nt = 0;
    long *ibase, *tibase, **linkage, **tlinkage;

    for (i = 1; i <= num_residue; i++)
        if (res_type[i] >= 0)
            num_nt++;  /* number of bases */
    if (!num_nt) {
        fprintf(stderr, "No base residues in this structure\n");
        return;
    }

    sprintf(filename, "%s%s", BDIR, "Block_R.alc");
    get_alc_nums(filename, &num, &nbond);

    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 1, 2);

    tnum = num * num_nt;
    tnbond = nbond * num_nt;

    tAtomName = cmatrix(1, tnum, 0, 2);
    txyz = dmatrix(1, tnum, 1, 3);
    tibase = lvector(1, tnum);
    tlinkage = lmatrix(1, tnbond, 1, 2);

    for (i = 1; i <= num_residue; i++) {
        if (res_type[i] == 1)  /* purine */
            sprintf(filename, "%s%s", BDIR, "Block_R.alc");
        else if (res_type[i] == 0)  /* pyrimidine */
            sprintf(filename, "%s%s", BDIR, "Block_Y.alc");
        else  /* not a nucleotide */
            continue;
        read_alc(filename, &num, &nbond, AtomName, xyz, ibase, linkage);

        base_idx(i, &bseq[i], &bidx, 1);  /* for coloring purpose */

        for (j = 1; j <= num; j++) {
            strcpy(tAtomName[ia + j], AtomName[j]);
            tibase[ia + j] = bidx;
            align_block_xyz(orien[i], org[i], xyz[j], txyz[ia + j]);
        }

        for (j = 1; j <= nbond; j++)
            for (k = 1; k <= 2; k++)
                tlinkage[ib + j][k] = ia + linkage[j][k];

        ia += num;
        ib += nbond;
    }
    write_alc(tnum, tnbond, tAtomName, txyz, tibase, tlinkage, alcfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_alc(tnum, tnbond, tAtomName, txyz, tibase, 1, tlinkage);
}

static void reset_alist_symbol_idx(char *alist, long *aidx)
{
    char c, asym[3], atoms_list[NELE][3];
    long i, k, nchar, natom;

    nchar = (long) strlen(alist);
    if (nchar % 2 != 0)
        fatal("each atomic symbol must be two-char long <%s>\n", alist);
    if (nchar == 0)
        fatal("empty atom list not allowed\n");

    upperstr(alist);  /* to all upper case */
    for (i = 0; i < nchar; i++) {
        c = alist[i];
        if ((c != '.') && !isupper((int) c))
            fatal("invalid char '%c': must be [.A-Z]\n", c);
        if (c == '.')
            alist[i] = ' ';
    }

    natom = nchar / 2;
    aidx[0] = natom;  /* number of atoms in the list */

    atom_info(1, atoms_list, NULL, NULL);

    for (i = 1; i <= natom; i++) {
        k = (i - 1) * 2;  /* 0, 2, 4 ... index position for next atom symbol */
        strncpy(asym, alist + k, 2);
        asym[2] = '\0';

        if (!num_strmatch(asym, Gvars.ATOM_NAMES, 0, Gvars.NUM_ELE))
            fatal("invalid atom symbol <%s>\n", asym);

        aidx[i] = asym_idx(asym, atoms_list, 0);
    }
}

static void check_space_in_altlist(char *alt_list)
{
    char str[BUF512];

    if (strchr(alt_list, ' ') == NULL) {  /* not included yet */
        sprintf(str, " %s", alt_list);
        strcpy(alt_list, str);
    }
}

static void check_if_valid_xml_parfile(void)
{
    static char *pars[] = {
        "min_base_hb",
        "hb_lower",
        "hb_dist1",
        "hb_dist2",
        "hb_atoms",
        "alt_list",
        "max_dorg",
        "min_dorg",
        "max_dv",
        "min_dv",
        "max_plane_angle",
        "min_plane_angle",
        "max_dNN",
        "min_dNN",
        "helix_break",
        "std_curved",
        "water_dist",
        "water_dlow",
        "water_atoms",
        "o3p_dist"
    };
    char *prefix = "", *dft = "xxxxxx", str[BUF512];
    char *p0, *line, fileLoc[BUF512];
    long i, num_pars = sizeof pars / sizeof pars[0];
    FILE *fp;

    get_BDIR(fileLoc, PAR_FILE);
    strcat(fileLoc, PAR_FILE);

    fp = open_file(fileLoc, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        line = trim(p0);  /* keep the original value of p0 */
        if (is_skip_line(line)) {
            free(p0);
            continue;
        }
        for (i = 0; i < num_pars; i++) {
            strcpy(str, dft);
            get_xml_tag(fp, prefix, line, " ", pars[i], str);
            if (!is_equal_string(str, dft)) {  /* with at least one valid xml tag */
                free(p0);
                close_file(fp);
                return;
            }
        }
        free(p0);
    }

    close_file(fp);

    fprintf(stderr, "File <%s> does NOT contain any valid xml tags for parameters\n"
            "\tmaybe it is in the previous v1.5 version?\n"
            "\tso system default setting is used here...\n", fileLoc);
}

static void mask_off_options(miscPars * misc_pars)
{
    misc_pars->min_base_hb = 1;
    misc_pars->hb_lower = 0.0;
    misc_pars->min_dorg = 0.0;
    misc_pars->min_dv = 0.0;
    misc_pars->min_plane_angle = 0.0;
    misc_pars->max_dNN = XBIG;
    misc_pars->water_dlow = 0.0;
}

void set_default_misc_pars(miscPars * misc_pars)
{
    misc_pars->min_base_hb = 1;  /* at least one H-bond between base atoms */
    /* 1jgq/1jgo/1jgp/1gix has N2 * N4  1.77 for pair C1_G-**+-C_C2 */
    misc_pars->hb_lower = 1.8;
    misc_pars->hb_dist1 = 4.0;
    misc_pars->hb_dist2 = 0.0;
    strcpy(misc_pars->hb_atoms, ".O.N");
    strcpy(misc_pars->alt_list, "A1");

    misc_pars->max_dorg = 15.0;
    misc_pars->min_dorg = 0.0;
    misc_pars->max_dv = 2.5;
    misc_pars->min_dv = 0.0;
    misc_pars->max_plane_angle = 65.0;
    misc_pars->min_plane_angle = 0.0;
    misc_pars->max_dNN = XBIG;
    misc_pars->min_dNN = 4.5;

    misc_pars->helix_break = 7.5;
    misc_pars->std_curved = 0.6;

    misc_pars->water_dist = 3.2;
    misc_pars->water_dlow = 0.0;
    strcpy(misc_pars->water_atoms, ".O.N");

    misc_pars->o3p_dist = 4.5;
}

static void output_misc_pars(FILE * fp, miscPars * misc_pars, char sep)
{
    long k = 76;

    print_sep(fp, sep, k);

    fprintf(fp, "<min_base_hb>%ld</min_base_hb>\n", misc_pars->min_base_hb);
    fprintf(fp, "<hb_lower>%g</hb_lower>\n", misc_pars->hb_lower);
    fprintf(fp, "<hb_dist1>%g</hb_dist1>\n", misc_pars->hb_dist1);
    fprintf(fp, "<hb_dist2>%g</hb_dist2>\n", misc_pars->hb_dist2);
    fprintf(fp, "<hb_atoms>%s</hb_atoms>\n", misc_pars->hb_atoms);
    fprintf(fp, "<alt_list>%s</alt_list>\n", misc_pars->alt_list);

    fprintf(fp, "<max_dorg>%g</max_dorg>\n", misc_pars->max_dorg);
    fprintf(fp, "<min_dorg>%g</min_dorg>\n", misc_pars->min_dorg);
    fprintf(fp, "<max_dv>%g</max_dv>\n", misc_pars->max_dv);
    fprintf(fp, "<min_dv>%g</min_dv>\n", misc_pars->min_dv);
    fprintf(fp, "<max_plane_angle>%g</max_plane_angle>\n", misc_pars->max_plane_angle);
    fprintf(fp, "<min_plane_angle>%g</min_plane_angle>\n", misc_pars->min_plane_angle);
    fprintf(fp, "<max_dNN>%g</max_dNN>\n", misc_pars->max_dNN);
    fprintf(fp, "<min_dNN>%g</min_dNN>\n", misc_pars->min_dNN);

    fprintf(fp, "<helix_break>%g</helix_break>\n", misc_pars->helix_break);
    fprintf(fp, "<std_curved>%g</std_curved>\n", misc_pars->std_curved);

    fprintf(fp, "<water_dist>%g</water_dist>\n", misc_pars->water_dist);
    fprintf(fp, "<water_dlow>%g</water_dlow>\n", misc_pars->water_dlow);
    fprintf(fp, "<water_atoms>%s</water_atoms>\n", misc_pars->water_atoms);

    fprintf(fp, "<o3p_dist>%g</o3p_dist>\n", misc_pars->o3p_dist);

    print_sep(fp, sep, k);
    fprintf(fp, "\n");
}

static void check_misc_pars(miscPars * misc_pars)
{
    if (misc_pars->min_base_hb < 0 || misc_pars->min_base_hb > 3) {
        fprintf(stderr, "min_base_hb: <%ld> -- reset to the default (1)\n",
                misc_pars->min_base_hb);
        misc_pars->min_base_hb = 1;
    }

    if (misc_pars->hb_lower < 0.0)
        fatal("hb_lower (%g) < 0.0\n", misc_pars->hb_lower);

    if (misc_pars->hb_dist1 < 0.0)
        fatal("hb_dist1 (%g) < 0.0\n", misc_pars->hb_dist1);

    if (misc_pars->hb_dist2 < 0.0)
        fatal("hb_dist2 (%g) < 0.0\n", misc_pars->hb_dist2);

    reset_alist_symbol_idx(misc_pars->hb_atoms, misc_pars->hb_idx);
    check_space_in_altlist(misc_pars->alt_list);

    if (misc_pars->max_dorg < 0.0)
        fatal("max_dorg (%g) < 0.0\n", misc_pars->max_dorg);

    if (misc_pars->min_dorg < 0.0)
        fatal("min_dorg (%g) < 0.0\n", misc_pars->min_dorg);
    if (misc_pars->min_dorg > misc_pars->max_dorg)
        fatal("min_dorg (%g) > max_dorg (%g)\n", misc_pars->min_dorg, misc_pars->max_dorg);

    if (misc_pars->max_dv < 0.0)
        fatal("max_dv (%g) < 0.0\n", misc_pars->max_dv);

    if (misc_pars->min_dv < 0.0)
        fatal("min_dv (%g) < 0.0\n", misc_pars->min_dv);
    if (misc_pars->min_dv > misc_pars->max_dv)
        fatal("min_dv (%g) > max_dv (%g)\n", misc_pars->min_dv, misc_pars->max_dv);

    if ((misc_pars->max_plane_angle < 0.0) || (misc_pars->max_plane_angle > 90.0))
        fatal("max_plane_angle (%g) must be in range [0, 90]\n", misc_pars->max_plane_angle);

    if ((misc_pars->min_plane_angle < 0.0) || (misc_pars->min_plane_angle > 90.0))
        fatal("min_plane_angle (%g) must be in range [0, 90]\n", misc_pars->min_plane_angle);
    if (misc_pars->min_plane_angle > misc_pars->max_plane_angle)
        fatal("min_plane_angle (%g) > max_plane_angle (%g)\n",
              misc_pars->min_plane_angle, misc_pars->max_plane_angle);

    if (misc_pars->max_dNN < 0.0)
        fatal("max_dNN (%g) < 0.0\n", misc_pars->max_dNN);

    if (misc_pars->min_dNN < 0.0)
        fatal("min_dNN (%g) < 0.0\n", misc_pars->min_dNN);
    if (misc_pars->min_dNN > misc_pars->max_dNN)
        fatal("min_dNN (%g) > max_dNN (%g)\n", misc_pars->min_dNN, misc_pars->max_dNN);

    if (misc_pars->helix_break < 0.0)
        fatal("helix_break (%g) < 0.0\n", misc_pars->helix_break);

    if (misc_pars->std_curved < 0.0)
        fatal("std_curved (%g) < 0.0\n", misc_pars->std_curved);

    if (misc_pars->water_dist < 0.0)
        fatal("water_dist (%g) < 0.0\n", misc_pars->water_dist);

    if (misc_pars->water_dlow < 0.0)
        fatal("water_dlow (%g) < 0.0\n", misc_pars->water_dlow);
    if (misc_pars->water_dlow > misc_pars->water_dist)
        fatal("water_dlow (%g) > water_dist (%g)\n", misc_pars->water_dlow, misc_pars->water_dist);

    reset_alist_symbol_idx(misc_pars->water_atoms, misc_pars->water_idx);

    if (misc_pars->o3p_dist < 0.0)
        fatal("o3p_dist (%g) < 0.0\n", misc_pars->o3p_dist);
}

static void get_3dna_pars(miscPars * misc_pars)
{
    char *p0, *line, *prefix = "", fileLoc[BUF512];
    FILE *fp, *fp2 = NULL;

    check_if_valid_xml_parfile();
    set_default_misc_pars(misc_pars);

    if (Gvars.DEBUG >= DEBUG_LEVEL) {
        fp2 = open_file("chk_misc_3dna.par", "w");
        output_misc_pars(fp2, misc_pars, '0');
    }

    get_BDIR(fileLoc, PAR_FILE);
    strcat(fileLoc, PAR_FILE);
    if (Gvars.VERBOSE)
        fprintf(stderr, " ...... reading file: %s ...... \n", PAR_FILE);

    fp = open_file(fileLoc, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        line = trim(p0);  /* keep the original value of p0 */
        if (is_skip_line(line)) {
            free(p0);
            continue;
        }

        get_xml_tag_long(fp, prefix, line, "min_base_hb", &misc_pars->min_base_hb);
        get_xml_tag_double(fp, prefix, line, "hb_lower", &misc_pars->hb_lower);
        get_xml_tag_double(fp, prefix, line, "hb_dist1", &misc_pars->hb_dist1);
        get_xml_tag_double(fp, prefix, line, "hb_dist2", &misc_pars->hb_dist2);
        get_xml_tag(fp, prefix, line, " ", "hb_atoms", misc_pars->hb_atoms);
        get_xml_tag(fp, prefix, line, " ", "alt_list", misc_pars->alt_list);

        get_xml_tag_double(fp, prefix, line, "max_dorg", &misc_pars->max_dorg);
        get_xml_tag_double(fp, prefix, line, "min_dorg", &misc_pars->min_dorg);
        get_xml_tag_double(fp, prefix, line, "max_dv", &misc_pars->max_dv);
        get_xml_tag_double(fp, prefix, line, "min_dv", &misc_pars->min_dv);
        get_xml_tag_double(fp, prefix, line, "max_plane_angle", &misc_pars->max_plane_angle);
        get_xml_tag_double(fp, prefix, line, "min_plane_angle", &misc_pars->min_plane_angle);
        get_xml_tag_double(fp, prefix, line, "max_dNN", &misc_pars->max_dNN);
        get_xml_tag_double(fp, prefix, line, "min_dNN", &misc_pars->min_dNN);

        get_xml_tag_double(fp, prefix, line, "helix_break", &misc_pars->helix_break);
        get_xml_tag_double(fp, prefix, line, "std_curved", &misc_pars->std_curved);

        get_xml_tag_double(fp, prefix, line, "water_dist", &misc_pars->water_dist);
        get_xml_tag_double(fp, prefix, line, "water_dlow", &misc_pars->water_dlow);
        get_xml_tag(fp, prefix, line, " ", "water_atoms", misc_pars->water_atoms);

        get_xml_tag_double(fp, prefix, line, "o3p_dist", &misc_pars->o3p_dist);

        free(p0);
    }

    close_file(fp);

    if (Gvars.DEBUG >= DEBUG_LEVEL)
        output_misc_pars(fp2, misc_pars, '1');

    check_misc_pars(misc_pars);

    if (Gvars.DEBUG >= DEBUG_LEVEL) {
        output_misc_pars(fp2, misc_pars, '2');
        close_file(fp2);
    }

    if (Gvars.DEBUG < DEBUG_LEVEL)
        mask_off_options(misc_pars);
}

static long overwrite_misc_pars(char *option)
{
    miscPars *misc_pars = &Gvars.misc_pars;

    if (case_str_pmatch(option, "-min_base_hb")) {
        misc_pars->min_base_hb = get_lvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-hb_lower")) {
        misc_pars->hb_lower = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-hb_dist1")) {
        misc_pars->hb_dist1 = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-hb_dist2")) {
        misc_pars->hb_dist2 = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-hb_atoms")) {
        get_strvalue(option, misc_pars->hb_atoms, FALSE);
        reset_alist_symbol_idx(misc_pars->hb_atoms, misc_pars->hb_idx);
        return TRUE;

    } else if (case_str_pmatch(option, "-alt_list")) {
        get_strvalue(option, misc_pars->alt_list, FALSE);
        check_space_in_altlist(misc_pars->alt_list);
        return TRUE;

    } else if (case_str_pmatch(option, "-max_dorg")) {
        misc_pars->max_dorg = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-min_dorg")) {
        misc_pars->min_dorg = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-max_dv")) {
        misc_pars->max_dv = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-min_dv")) {
        misc_pars->min_dv = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-max_plane_angle")) {
        misc_pars->max_plane_angle = get_dvalue(option, 0, 90);
        return TRUE;

    } else if (case_str_pmatch(option, "-min_plane_angle")) {
        misc_pars->min_plane_angle = get_dvalue(option, 0, 90);
        return TRUE;

    } else if (case_str_pmatch(option, "-max_dNN")) {
        misc_pars->max_dNN = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-min_dNN")) {
        misc_pars->min_dNN = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-helix_break")) {
        misc_pars->helix_break = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-std_curved")) {
        misc_pars->std_curved = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-water_dist")) {
        misc_pars->water_dist = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-water_dlow")) {
        misc_pars->water_dlow = get_dvalue(option, 0, 999);
        return TRUE;

    } else if (case_str_pmatch(option, "-water_atoms")) {
        get_strvalue(option, misc_pars->water_atoms, FALSE);
        reset_alist_symbol_idx(misc_pars->water_atoms, misc_pars->water_idx);
        return TRUE;

    } else if (case_str_pmatch(option, "-o3p_dist")) {
        misc_pars->o3p_dist = get_dvalue(option, 0, 999);
        return TRUE;
    }

    return FALSE;
}

void lsplane_xyz(double **xyz, long num_plane, long *atom_plane, double **nxyz, double *z)
{
    long i, k;
    double **xyz_plane;
    double odist, adist[BUF512], ppos[4];

    xyz_plane = dmatrix(1, num_plane, 1, 3);

    for (i = 1; i <= num_plane; i++) {
        k = atom_plane[i];
        cpxyz(xyz[k], xyz_plane[i]);
    }

    ls_plane(xyz_plane, num_plane, z, ppos, &odist, adist);
    plane_xyz(num_plane, xyz_plane, ppos, z, nxyz);

    free_dmatrix(xyz_plane, 1, num_plane, 1, 3);
}

/* base-pairing type information as in check_pair():
 * +1 WC geometry; +2: WC pair; -1: other cases */
long read_PairInfo(char *inpfile, long **pair_info)
{
    char *p0, *p;
    long i, num_bp;
    FILE *fp;

    fp = open_file(inpfile, "r");

    for (i = 1; i <= 5; i++) {  /* skip top 5 lines */
        p0 = my_getline(fp);
        if (i == 4 && sscanf(p0, "%ld", &num_bp) != 1)
            fatal("\tcan't extract number of base-pairs <%s>\n", p0);
        free(p0);
    }

    for (i = 1; i <= num_bp; i++) {
        p0 = my_getline(fp);
        if (sscanf(p0, "%ld %ld", &pair_info[i][1], &pair_info[i][2]) != 2)
            fatal("\tcan't extract pairing number %ld <%s>\n", i, p0);
        pair_info[i][0] = -1;
        p = strchr(p0, ']') + 3;  /* ]G-----C[ */
        if (*(p + 1) == '-')  /* WC geometry */
            pair_info[i][0] = 1;
        if (*p == '-')  /* WC pair */
            pair_info[i][0] = 2;
        free(p0);
    }

    close_file(fp);

    return num_bp;
}

/* use O3' and P distance to check if residues i & j are connected by gap-residues */
long is_linked_by_gap(long i, long j, double **o3_p)
{
    long k;

    for (k = i; k < j; k++)
        if (!is_linked(k, k + 1, o3_p))
            return 0L;

    return 1L;
}

/* use O3' and P distance to check if residues i & j are directly connected */
long is_linked(long i, long j, double **o3_p)
{
    double d;

    d = distance_ab(o3_p, i, j, 4, 8);  /* O3'[i]-P[j] */
    if (dval_in_range(d, 0.0, O3P_UPPER))
        return 1L;  /* 5'--->3' linkage */

    d = distance_ab(o3_p, j, i, 4, 8);  /* O3'[j]-P[i] */
    if (dval_in_range(d, 0.0, O3P_UPPER))
        return -1L;  /* 3'--->5' linkage */

    return 0L;  /* no direct linkage */
}

/* calculate distance between ia & ib: ipa & ipb mark if ia/ib exist */
double distance_ab(double **o3_p, long ia, long ib, long ipa, long ipb)
{
    return (o3_p[ia][ipa] > 0.0 && o3_p[ib][ipb] > 0.0) ?
        p1p2_dist(o3_p[ia] + ipa - 4, o3_p[ib] + ipb - 4) : -1.0;
}

/* write the position and orientation of the transformation to ROTMAT_FILE
 *     To be used with rotate_mol with -l and -t=ROTMAT option */
void write_rotmat(double **rotmat)
{
    char *format = "%12.4f%12.4f%12.4f\n";
    long i;
    FILE *fp;

    fp = open_file(ROTMAT_FILE, "w");
    fprintf(fp, "    1  # x-, y-, z-axes row-rise\n");
    fprintf(fp, format, rotmat[4][1], rotmat[4][2], rotmat[4][3]);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, rotmat[1][i], rotmat[2][i], rotmat[3][i]);

    close_file(fp);
}

/* read in a rotation matrix: x-axis can be any vector;
   y-axis orthogonal to it; z-axis forms a right-handed system */
void read_rotmat(char *rotfile, double **rotmat)
{
    char str[BUF512], *format = "%lf %lf %lf";
    double tmp[4], **tmprot;
    long i, row_wise;
    FILE *fp;

    fp = open_file(rotfile, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", &row_wise) != 1)
        fatal("error reading rotmat file: row/column-wise indicator\n");
    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, format, &rotmat[4][1], &rotmat[4][2], &rotmat[4][3]) != 3)
        fatal("error reading origin line\n");

    tmprot = dmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++) {
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &tmp[1], &tmp[2], &tmp[3]) != 3)
            fatal("error reading x/y/z directions\n");
        vec_norm(tmp);  /* normalize it */
        cpxyz(tmp, tmprot[i]);
    }
    close_file(fp);

    if (!row_wise) {  /* make tmprot row-vector for xyz axes */
        transpose_matrix(tmprot, 3, 3, rotmat);
        copy_dmatrix(rotmat, 3, 3, tmprot);
    }
    if (fabs(dot(tmprot[1], tmprot[2])) >= 0.005)
        fprintf(stderr, "Warning: x- and y-axes are not orthogonal\n");
    vec_orth(tmprot[2], tmprot[1]);  /* orthogonal component to x-axis as y-axis */
    cross(tmprot[1], tmprot[2], tmp);
    if (fabs(dot(tmp, tmprot[3])) <= 0.9995)
        fprintf(stderr, "Warning: z-axis is not perpendicular to xy-plane\n");
    cpxyz(tmp, tmprot[3]);
    transpose_matrix(tmprot, 3, 3, rotmat);  /* make rotmat column-wise */

    free_dmatrix(tmprot, 1, 3, 1, 3);
}

void reverse_y_z_columns(double **R)
{
    long i;

    for (i = 1; i <= 3; i++) {
        R[i][2] = -R[i][2];  /* reverse y-axis */
        R[i][3] = -R[i][3];  /* reverse z-axis */
    }
}

long get_num_nt(long num_residue, long *RY)
{
    long i, num_nt = 0;

    for (i = 1; i <= num_residue; i++)
        if (RY[i] >= 0)
            num_nt++;

    return num_nt;
}

/* get the local reference frame for each peptide unit defined by 6 atoms:
   CA(i), C(i), O(i), N(i+1), H(i+1) & CA(i+1) -- NB: Hs normally N/A */
void peptide_frame(long num_residue, char *BDIR, long *res_type, long *mchain,
                   double **xyz, double **orien, double **org)
{
    char filename[BUF512], *sChainID, **sAtomName, **sResName;
    double rms_fit, **ePxyz, **fitted_xyz, **sPxyz, **sxyz, **R;
    long i, ioffset, j, nmatch, pnum = 6;  /* 6 atoms in snap_pep.pdb */
    long idx[7], *sResSeq;

    sprintf(filename, "%s%s", BDIR, SNAP_PEP_PDB);
    fprintf(stderr, " ...... reading file: %s ...... \n", filename);

    sAtomName = cmatrix(1, pnum, 0, 4);  /* for standard geometry */
    sResName = cmatrix(1, pnum, 0, 3);
    sChainID = cvector(1, pnum);
    sResSeq = lvector(1, pnum);
    sxyz = dmatrix(1, pnum, 1, 3);
    read_pdb(filename, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL, 1, "*");

    ePxyz = dmatrix(1, pnum, 1, 3);
    sPxyz = dmatrix(1, pnum, 1, 3);
    fitted_xyz = dmatrix(1, pnum, 1, 3);
    R = dmatrix(1, 3, 1, 3);

    org[num_residue][4] = -1.0;
    for (i = 1; i < num_residue; i++) {  /* between i and i + 1 */
        org[i][4] = -1.0;  /* as a default indicator */
        if (res_type[i] != -1 || res_type[i + 1] != -1)  /* non-amino-acid */
            continue;
        ioffset = (i - 1) * 4;  /* for residue i */
        idx[1] = mchain[ioffset + 2];  /* CA (i) */
        idx[2] = mchain[ioffset + 3];  /* C (i) */
        idx[3] = mchain[ioffset + 4];  /* O (i) */
        ioffset += 4;  /* for residue i + 1 */
        idx[4] = mchain[ioffset + 1];  /* N (i + 1) */
        idx[5] = 0;  /* H (i + 1): not counted */
        idx[6] = mchain[ioffset + 2];  /* CA (i + 1) */

        if (!idx[2] || !idx[4])  /* C(i) or N(i + 1) does not exist */
            continue;
        if (p1p2_dist(xyz[idx[2]], xyz[idx[4]]) > 2.0)  /* 1.32 C-N */
            continue;  /* AAs (i) and (i + 1) NOT connected by a peptide bond */

        nmatch = 0;
        for (j = 1; j <= pnum; j++) {
            if (idx[j]) {
                ++nmatch;
                cpxyz(xyz[idx[j]], ePxyz[nmatch]);
                cpxyz(sxyz[j], sPxyz[nmatch]);
            }
        }

        rms_fit = ls_fitting(sPxyz, ePxyz, nmatch, fitted_xyz, R, org[i]);

        org[i][4] = rms_fit;
        mst2orien(orien[i], 0, R);  /* attached to residue i */

        if (Gvars.VERBOSE)
            fprintf(stderr, "RMS for peptide plane %ld: %g\n", i, rms_fit);
    }

    free_pdb(pnum, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL);
    free_dmatrix(ePxyz, 1, pnum, 1, 3);
    free_dmatrix(sPxyz, 1, pnum, 1, 3);
    free_dmatrix(fitted_xyz, 1, pnum, 1, 3);
    free_dmatrix(R, 1, 3, 1, 3);
}

/* representing each peptide bond as a rectangular block */
void peptide_blks(long num_residue, char *BDIR, long *cidx, double **orien, double **org,
                  char *alcfile)
{
    char filename[BUF512], **AtomName, **tAtomName;
    double **txyz, **xyz;
    long i, j, k, nbond, num, tnbond, tnum;
    long ia = 0, ib = 0, num_pep = 0;
    long *ibase, *tibase, **linkage, **tlinkage;

    for (i = 1; i < num_residue; i++)  /* NB: org[i][4] default to -1.0 */
        if (dval_in_range(org[i][4], 0.0, TRSP_RMS))
            num_pep++;
    if (!num_pep) {
        fprintf(stderr, "No trans peptide bonds in this structure\n");
        return;
    }

    sprintf(filename, "%s%s", BDIR, SNAP_PEP_ALC);
    fprintf(stderr, " ...... reading file: %s ...... \n", filename);

    get_alc_nums(filename, &num, &nbond);
    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);  /* of no use here */
    linkage = lmatrix(1, nbond, 1, 2);
    read_alc(filename, &num, &nbond, AtomName, xyz, ibase, linkage);

    tnum = num * num_pep;  /* total # of vertices & linkages */
    tnbond = nbond * num_pep;

    tAtomName = cmatrix(1, tnum, 0, 2);
    txyz = dmatrix(1, tnum, 1, 3);
    tibase = lvector(1, tnum);
    tlinkage = lmatrix(1, tnbond, 1, 2);

    for (i = 1; i < num_residue; i++) {
        if (!dval_in_range(org[i][4], 0.0, TRSP_RMS))
            continue;  /* not a trans peptide bond */
        for (j = 1; j <= num; j++) {  /* for current block vertices */
            strcpy(tAtomName[ia + j], AtomName[j]);
            tibase[ia + j] = cidx[i];  /* cidx is for [1 .. num_residue] */
            align_block_xyz(orien[i], org[i], xyz[j], txyz[ia + j]);
        }
        for (j = 1; j <= nbond; j++)
            for (k = 1; k <= 2; k++)
                tlinkage[ib + j][k] = ia + linkage[j][k];
        ia += num;
        ib += nbond;
    }
    write_alc(tnum, tnbond, tAtomName, txyz, tibase, tlinkage, alcfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_alc(tnum, tnbond, tAtomName, txyz, tibase, 1, tlinkage);
}

/* this function must be run before calling atom_idx() */
void set_my_globals(char *pgname)
{

    Gvars.VERBOSE = FALSE;
    Gvars.CHAIN_CASE = TRUE;  /* case-sensitive */
    Gvars.ALL_MODEL = FALSE;  /* just one model: ENDMDL/END */
    Gvars.DEBUG = FALSE;  /* for general distribution */
    Gvars.PROGNAME = pgname;
    Gvars.ATTACH_RESIDUE = TRUE;  /* add connected HETATM metals or residues */
    Gvars.THREE_LETTER_NTS = FALSE;  /* ADE/CYT/GUA/THY/URA */
    Gvars.PDBV3 = TRUE;
    Gvars.ORIGINAL_COORDINATE = FALSE;
    Gvars.OCCUPANCY = FALSE;
    Gvars.HEADER = TRUE;

    Gvars.NT_CUTOFF = 0.2618;  /* 2o8b C30 on F has 0.24 */

    get_3dna_homedir(Gvars.X3DNA_HOMEDIR);
    get_3dna_version(Gvars.X3DNA_HOMEDIR, Gvars.X3DNA_VER);
    strcpy(Gvars.CHAIN_MARKERS, CMARKERS);
    strcpy(Gvars.REBUILD_CHAIN_IDS, REBUILD_CIDS);

    Gvars.ATOM_NAMES = ptable;
    Gvars.NUM_ELE = sizeof ptable / sizeof ptable[0] - 1;  /* minus - 1 */

    Gvars.ATOMLIST = cmatrix(1, BUFBIG, 0, 6);
    get_atomlist(Gvars.ATOMLIST, &Gvars.NUM_SATOM);

    Gvars.BASELIST = cmatrix(1, BUFBIG, 0, 4);
    get_baselist(Gvars.BASELIST, &Gvars.NUM_SBASE);

    Gvars.AtomName0 = NULL;
    Gvars.ResName0 = NULL;
    Gvars.Name0 = FALSE;

    get_3dna_pars(&Gvars.misc_pars);
}

void clear_my_globals(void)
{
    free_cmatrix(Gvars.ATOMLIST, 1, -1, 0, -1);
    free_cmatrix(Gvars.BASELIST, 1, -1, 0, -1);
}

static void set_chain_markers(char *option)
{
    long k;

    get_strvalue(option, Gvars.CHAIN_MARKERS, FALSE);
    k = strlen(Gvars.CHAIN_MARKERS);
    if (k < 4) {
        fprintf(stderr, "too short [< 4-char] input for option: %s\n", option);
        fprintf(stderr, "\treset to default: -chain_markers='%s'\n", CMARKERS);
        strcpy(Gvars.CHAIN_MARKERS, CMARKERS);
    } else if (k == 4) {
        if (Gvars.CHAIN_MARKERS[0] == CMARKERS[0])
            Gvars.CHAIN_MARKERS[4] = CMARKERS[4];
        else
            Gvars.CHAIN_MARKERS[4] = '1';
    }
}

static void set_rebuild_cids(char *option)
{
    char cids[BUF512];
    long k;

    strcpy(cids, "");
    get_strvalue(option, cids, FALSE);
    k = strlen(cids);

    if (k == 0)
        strcpy(cids, REBUILD_CIDS);
    else if (k == 1)
        strcat(cids, "B");

    cids[2] = '\0';  /* keep only 2 chars */
    fprintf(stderr, "[i] -- chain ids for rebuild: '%s'\n", cids);

    strcpy(Gvars.REBUILD_CHAIN_IDS, cids);
}

long check_global_options(char *option)
{
    if (lux_ncmatch(option, "^--?debug")) {
        Gvars.DEBUG = get_lvalue(option, 0, BUF512);
        return TRUE;

    } else if (lux_ncmatch(option, "^--?chain_case")) {
        Gvars.CHAIN_CASE = set_switch_default_true(option);
        return TRUE;

    } else if (lux_ncmatch(option, "^--?chain_markers")) {
        set_chain_markers(option);
        return TRUE;

    } else if (lux_ncmatch(option, "^--?rebuild-cids")) {
        set_rebuild_cids(option);
        return TRUE;

    } else if (lux_ncmatch(option, "^--?(all_model|models|symm)")) {
        Gvars.ALL_MODEL = set_switch_default_true(option);
        return TRUE;

    } else if (lux_ncmatch(option, "^--?verbose")) {
        Gvars.VERBOSE = set_switch_default_true(option);
        return TRUE;

    } else if (lux_ncmatch(option, "^--?attach")) {
        Gvars.ATTACH_RESIDUE = set_switch_default_true(option);
        return TRUE;

    } else if (lux_ncmatch(option, "^--?three")) {
        Gvars.THREE_LETTER_NTS = set_switch_default_true(option);
        return TRUE;

    } else if (lux_ncmatch(option, "^--?pdbv3")) {
        Gvars.PDBV3 = set_switch_default_true(option);
        return TRUE;

    } else if (lux_ncmatch(option, "^--?(ori|raw)")) {
        Gvars.ORIGINAL_COORDINATE = set_switch_default_true(option);
        return TRUE;

    } else if (lux_ncmatch(option, "^--?occ")) {
        Gvars.OCCUPANCY = set_switch_default_true(option);
        return TRUE;

    } else if (lux_ncmatch(option, "^--?hea")) {
        Gvars.HEADER = set_switch_default_true(option);
        if (Gvars.HEADER && lux_ncmatch(option, "^--?(whole|all"))
            Gvars.HEADER = BUF32;
        return TRUE;

    } else if (lux_ncmatch(option, "^--?ntc")) {
        Gvars.NT_CUTOFF = get_dvalue(option, 0.05, 1.0);
        return TRUE;

    } else if (overwrite_misc_pars(option))
        return TRUE;

    else
        return FALSE;
}

/* ++++++++++++++++++++++++++++++++++++ SNAP-related functions ++++++++++++++++++++++++++++++++++++ */

/* get the local AA reference frame based on Ca-Cb-N as in Pabo/Honig */
void get_AA_frames(long num_residue, long **seidx, long *res_type, char **AtomName,
                   char **ResName, char *ChainID, long *ResSeq, char **Miscs,
                   double **xyz, char *BDIR, double **orien, double **org)
{
    char idmsg[BUF512];
    char filename[BUF512], *sChainID, **sAtomName, **sResName;
    double rms_fit, **ePxyz, **fitted_xyz, **sPxyz, **sxyz, **R;
    long i, ib, ie, j, nmatch, pnum = 4;  /* 4 atoms in SNAP_AAA */
    long idx[5], *sResSeq;

    sprintf(filename, "%s%s", BDIR, SNAP_AAA);
    fprintf(stderr, " ...... reading file: %s ...... \n", filename);

    sAtomName = cmatrix(1, pnum, 0, 4);  /* for standard geometry */
    sResName = cmatrix(1, pnum, 0, 3);
    sChainID = cvector(1, pnum);
    sResSeq = lvector(1, pnum);
    sxyz = dmatrix(1, pnum, 1, 3);
    read_pdb(filename, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL, 1, "*");

    ePxyz = dmatrix(1, pnum, 1, 3);
    sPxyz = dmatrix(1, pnum, 1, 3);
    fitted_xyz = dmatrix(1, pnum, 1, 3);
    R = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= num_residue; i++) {
        if (res_type[i] != -1)  /* non-amino-acid */
            continue;
        ib = seidx[i][1];
        ie = seidx[i][2];
        get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);
        idx[1] = find_1st_atom(" N  ", AtomName, ib, ie, idmsg);
        idx[2] = find_1st_atom(" CA ", AtomName, ib, ie, idmsg);
        idx[3] = find_1st_atom(" C  ", AtomName, ib, ie, idmsg);
        idx[4] = find_1st_atom(" CB ", AtomName, ib, ie, "");  /* for glycine */

        org[i][4] = -1.0;  /* as a default indicator */
        nmatch = 0;
        for (j = 1; j <= pnum; j++) {
            if (idx[j]) {
                ++nmatch;
                cpxyz(xyz[idx[j]], ePxyz[nmatch]);
                cpxyz(sxyz[j], sPxyz[nmatch]);
            }
        }

        if (nmatch >= 3) {  /* extract match Cb xyz later ... */
            rms_fit = ls_fitting(sPxyz, ePxyz, nmatch, fitted_xyz, R, org[i]);
            org[i][4] = rms_fit;
            mst2orien(orien[i], 0, R);
            if (Gvars.VERBOSE) {
                fprintf(stderr, "RMS %s %g\n", idmsg, rms_fit);
                verify_Cb_coordinates(nmatch, ePxyz, fitted_xyz);
            }
        }
    }

    free_pdb(pnum, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL);
    free_dmatrix(ePxyz, 1, pnum, 1, 3);
    free_dmatrix(sPxyz, 1, pnum, 1, 3);
    free_dmatrix(fitted_xyz, 1, pnum, 1, 3);
    free_dmatrix(R, 1, 3, 1, 3);
}

void verify_Cb_coordinates(long nmatch, double **ePxyz, double **fitted_xyz)
{
    char *fmt = "\t%9.3f%9.3f%9.3f\n";

    if (nmatch == 4) {
        fprintf(stderr, fmt, ePxyz[4][1], ePxyz[4][2], ePxyz[4][3]);
        fprintf(stderr, fmt, fitted_xyz[4][1], fitted_xyz[4][2], fitted_xyz[4][3]);
    }
}

void residue_chain_resnum(char chain_id, long res_seq, char *misc, char *idmsg)
{
    char iCode = misc[2];

    sprintf(idmsg, "%c%ld%c", chain_id, res_seq, iCode);
}

void residue_strid(char chain_id, long res_seq, char *misc, char *rname, char *idmsg)
{
    char modelNum[10], iCode = misc[2], newName[BUF512];
    long model_lval, nlen = 4;

    strncpy(modelNum, misc + 30, nlen);
    modelNum[nlen] = '\0';
    model_lval = cvt2long(modelNum);  /* model number of the residue */

    if (iCode == ' ')
        iCode = '.';
    convert_resNameSpace(rname, '_', newName);  /* '\0' to eliminate space */

    if (model_lval == 0)
        sprintf(idmsg, "%c%ld%c%s", chain_id, res_seq, iCode, newName);
    else
        sprintf(idmsg, "%4s%c%ld%c%s", modelNum, chain_id, res_seq, iCode, newName);
}

void convert_resNameSpace(char *resName, char replacement, char *newName)
{
    long i, j = 0, k = strlen(resName);

    if (replacement != '\0') {
        for (i = 0; i < k; i++)
            newName[i] = (resName[i] == ' ') ? replacement : resName[i];
        newName[i] = '\0';
    } else {  /* eliminate ' ' */
        for (i = 0; i < k; i++)
            if (resName[i] != ' ')
                newName[j++] = resName[i];
        newName[j] = '\0';
    }
}

void get_planarFrame(char *aa, long pnum, char *planar_atoms[], char *BDIR, char *idmsg,
                     long ib, long ie, char **AtomName, double **xyz, double *orien_i,
                     double *org_i)
/* see also function set_planarStr() in set_planar.c */
{
    char spdb[BUF512], sidmsg[BUF512], *sChainID, **sAtomName, **sResName;
    long i, exp_katom, std_katom, nmatch = 0, snum, *sResSeq;
    double rms_fit, **ePxyz, **fitted_xyz, **sPxyz, **sxyz, **R;

    sprintf(spdb, "%s%s%s%s", BDIR, "snap_", aa, ".pdb");
    sprintf(sidmsg, "in standard residue file %s", spdb);

    sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
    sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
    sChainID = cvector(1, NUM_RESIDUE_ATOMS);
    sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
    sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);

    snum = read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL, 1, "*");

    ePxyz = dmatrix(1, pnum, 1, 3);
    sPxyz = dmatrix(1, pnum, 1, 3);
    fitted_xyz = dmatrix(1, pnum, 1, 3);
    R = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= pnum; i++) {
        exp_katom = find_1st_atom(planar_atoms[i], AtomName, ib, ie, idmsg);
        std_katom = find_1st_atom(planar_atoms[i], sAtomName, 1, snum, sidmsg);
        if (exp_katom && std_katom) {
            ++nmatch;
            cpxyz(xyz[exp_katom], ePxyz[nmatch]);
            cpxyz(sxyz[std_katom], sPxyz[nmatch]);
        }
    }

    if (nmatch >= 3) {
        rms_fit = ls_fitting(sPxyz, ePxyz, nmatch, fitted_xyz, R, org_i);
        org_i[4] = rms_fit;
        mst2orien(orien_i, 0, R);
        if (Gvars.VERBOSE)
            fprintf(stderr, "%s ===> %g\n", idmsg, org_i[4]);
    }

    free_pdb(NUM_RESIDUE_ATOMS, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL);
    free_dmatrix(ePxyz, 1, pnum, 1, 3);
    free_dmatrix(sPxyz, 1, pnum, 1, 3);
    free_dmatrix(fitted_xyz, 1, pnum, 1, 3);
    free_dmatrix(R, 1, 3, 1, 3);
}

void get_argFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i)
{
    static char *planar_atoms[] = { "XXXX", " CD ", " NE ", " CZ ", " NH1", " NH2" };
    long pnum = sizeof planar_atoms / sizeof planar_atoms[0] - 1;

    get_planarFrame("arg", pnum, planar_atoms, BDIR, idmsg, ib, ie, AtomName, xyz, orien_i, org_i);
}

void get_pheFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i)
{
    static char *planar_atoms[] =
        { "XXXX", " CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ " };
    long pnum = sizeof planar_atoms / sizeof planar_atoms[0] - 1;

    get_planarFrame("phe", pnum, planar_atoms, BDIR, idmsg, ib, ie, AtomName, xyz, orien_i, org_i);
}

void get_tyrFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i)
{
    static char *planar_atoms[] =
        { "XXXX", " CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ ", " OH " };
    long pnum = sizeof planar_atoms / sizeof planar_atoms[0] - 1;

    get_planarFrame("tyr", pnum, planar_atoms, BDIR, idmsg, ib, ie, AtomName, xyz, orien_i, org_i);
}

void get_trpFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i)
{
    static char *planar_atoms[] = { "XXXX", " CB ", " CG ", " CD1", " CD2", " NE1", " CE2", " CE3",
        " CZ2", " CZ3", " CH2"
    };
    long pnum = sizeof planar_atoms / sizeof planar_atoms[0] - 1;

    get_planarFrame("trp", pnum, planar_atoms, BDIR, idmsg, ib, ie, AtomName, xyz, orien_i, org_i);
}

void get_hisFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i)
{
    static char *planar_atoms[] = { "XXXX", " CB ", " CG ", " ND1", " CD2", " CE1", " NE2" };
    long pnum = sizeof planar_atoms / sizeof planar_atoms[0] - 1;

    get_planarFrame("his", pnum, planar_atoms, BDIR, idmsg, ib, ie, AtomName, xyz, orien_i, org_i);
}

void get_lysFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i)
{
    static char *planar_atoms[] = { "XXXX", " CD ", " CE ", " NZ " };
    long pnum = sizeof planar_atoms / sizeof planar_atoms[0] - 1;

    if (!Gvars.VERBOSE)
        strcpy(idmsg, "");  /* no message of missing atoms for lysine */
    get_planarFrame("lys", pnum, planar_atoms, BDIR, idmsg, ib, ie, AtomName, xyz, orien_i, org_i);
}

void get_asnFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i)
{
    static char *planar_atoms[] = { "XXXX", " CB ", " CG ", " OD1", " ND2" };
    long pnum = sizeof planar_atoms / sizeof planar_atoms[0] - 1;

    get_planarFrame("asn", pnum, planar_atoms, BDIR, idmsg, ib, ie, AtomName, xyz, orien_i, org_i);
}

void get_glnFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i)
{
    static char *planar_atoms[] = { "XXXX", " CG ", " CD ", " OE1", " NE2" };
    long pnum = sizeof planar_atoms / sizeof planar_atoms[0] - 1;

    get_planarFrame("gln", pnum, planar_atoms, BDIR, idmsg, ib, ie, AtomName, xyz, orien_i, org_i);
}

void get_aspFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i)
{
    static char *planar_atoms[] = { "XXXX", " CB ", " CG ", " OD1", " OD2" };
    long pnum = sizeof planar_atoms / sizeof planar_atoms[0] - 1;

    get_planarFrame("asp", pnum, planar_atoms, BDIR, idmsg, ib, ie, AtomName, xyz, orien_i, org_i);
}

void get_gluFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i)
{
    static char *planar_atoms[] = { "XXXX", " CG ", " CD ", " OE1", " OE2" };
    long pnum = sizeof planar_atoms / sizeof planar_atoms[0] - 1;

    get_planarFrame("glu", pnum, planar_atoms, BDIR, idmsg, ib, ie, AtomName, xyz, orien_i, org_i);
}

void get_planarAA_frames(long num_residue, long **seidx, char **AtomName, char **ResName,
                         char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                         long *res_type, long *cidx)
{
    char idmsg[BUF512], BDIR[BUF512];
    long i, ib, ie;
    double **paa_orien, **paa_org;

    get_BDIR(BDIR, SNAP_PEP_PDB);

    paa_orien = dmatrix(1, num_residue, 1, 9);
    paa_org = dmatrix(1, num_residue, 1, 4);

    for (i = 1; i <= num_residue; i++) {
        paa_org[i][4] = -1;  /* as a default indicator */
        if (res_type[i] != -1)  /* non-amino-acid */
            continue;
        ib = seidx[i][1];
        ie = seidx[i][2];
        get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);
        if (strcmp(ResName[ib], "ARG") == 0)
            get_argFrame(BDIR, idmsg, ib, ie, AtomName, xyz, paa_orien[i], paa_org[i]);
        else if (strcmp(ResName[ib], "PHE") == 0)
            get_pheFrame(BDIR, idmsg, ib, ie, AtomName, xyz, paa_orien[i], paa_org[i]);
        else if (strcmp(ResName[ib], "TYR") == 0)
            get_tyrFrame(BDIR, idmsg, ib, ie, AtomName, xyz, paa_orien[i], paa_org[i]);
        else if (strcmp(ResName[ib], "TRP") == 0)
            get_trpFrame(BDIR, idmsg, ib, ie, AtomName, xyz, paa_orien[i], paa_org[i]);
        else if (strcmp(ResName[ib], "HIS") == 0)
            get_hisFrame(BDIR, idmsg, ib, ie, AtomName, xyz, paa_orien[i], paa_org[i]);
        else if (strcmp(ResName[ib], "LYS") == 0)
            get_lysFrame(BDIR, idmsg, ib, ie, AtomName, xyz, paa_orien[i], paa_org[i]);
        else if (strcmp(ResName[ib], "ASN") == 0)
            get_asnFrame(BDIR, idmsg, ib, ie, AtomName, xyz, paa_orien[i], paa_org[i]);
        else if (strcmp(ResName[ib], "GLN") == 0)
            get_glnFrame(BDIR, idmsg, ib, ie, AtomName, xyz, paa_orien[i], paa_org[i]);
        else if (strcmp(ResName[ib], "ASP") == 0)
            get_aspFrame(BDIR, idmsg, ib, ie, AtomName, xyz, paa_orien[i], paa_org[i]);
        else if (strcmp(ResName[ib], "GLU") == 0)
            get_gluFrame(BDIR, idmsg, ib, ie, AtomName, xyz, paa_orien[i], paa_org[i]);
    }

    planarAA_blks(num_residue, seidx, ResName, BDIR, cidx, paa_orien, paa_org, ABLKALC_FILE);

    free_dmatrix(paa_orien, 1, num_residue, 1, 9);
    free_dmatrix(paa_org, 1, num_residue, 1, 4);
}

/* representing the planar moiety of AA as a rectangular block */
void planarAA_blks(long num_residue, long **seidx, char **ResName, char *BDIR, long *cidx,
                   double **paa_orien, double **paa_org, char *alcfile)
{
    char aa[4], filename[BUF512], **AtomName, **tAtomName;
    double **txyz, **xyz;
    long i, j, k, nbond, num, tnbond, tnum;
    long ia = 0, ib = 0, num_paa = 0;
    long *ibase, *tibase, **linkage, **tlinkage;

    for (i = 1; i <= num_residue; i++)  /* NB: paa_org[i][4] default to -1.0 */
        if (dval_in_range(paa_org[i][4], 0.0, TRSP_RMS))
            num_paa++;
    if (!num_paa) {
        fprintf(stderr, "no planar amino acids in this structure\n");
        return;
    }

    sprintf(filename, "%s%s%s%s", BDIR, "snap_", "arg", ".alc");
    get_alc_nums(filename, &num, &nbond);

    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);  /* of no use here */
    linkage = lmatrix(1, nbond, 1, 2);

    tnum = num * num_paa;  /* total # of vertices & linkages */
    tnbond = nbond * num_paa;

    tAtomName = cmatrix(1, tnum, 0, 2);
    txyz = dmatrix(1, tnum, 1, 3);
    tibase = lvector(1, tnum);
    tlinkage = lmatrix(1, tnbond, 1, 2);

    for (i = 1; i <= num_residue; i++) {
        if (!dval_in_range(paa_org[i][4], 0.0, TRSP_RMS))
            continue;  /* not AA with planar moiety */
        strcpy(aa, ResName[seidx[i][1]]);
        lowerstr(aa);
        sprintf(filename, "%s%s%s%s", BDIR, "snap_", aa, ".alc");
        read_alc(filename, &num, &nbond, AtomName, xyz, ibase, linkage);

        for (j = 1; j <= num; j++) {  /* for current block vertices */
            strcpy(tAtomName[ia + j], AtomName[j]);
            tibase[ia + j] = cidx[i];  /* cidx is for [1 .. num_residue] */
            align_block_xyz(paa_orien[i], paa_org[i], xyz[j], txyz[ia + j]);
        }
        for (j = 1; j <= nbond; j++)
            for (k = 1; k <= 2; k++)
                tlinkage[ib + j][k] = ia + linkage[j][k];
        ia += num;
        ib += nbond;
    }
    write_alc(tnum, tnbond, tAtomName, txyz, tibase, tlinkage, alcfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_alc(tnum, tnbond, tAtomName, txyz, tibase, 1, tlinkage);
}

void print_resid(long num_residue, long **seidx, char **ResName, char *ChainID,
                 long *ResSeq, char **Miscs, long *res_type)
{
    long i, ib;

    if (!Gvars.VERBOSE)
        return;

    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        fprintf(stderr, "%5ld %2ld %c %s %4ld%c\n", i, res_type[i], ChainID[ib],
                ResName[ib], ResSeq[ib], Miscs[ib][2]);
    }
}

void snap_atype(char **AtomName, long num_residue, long **seidx, long *res_type, long **atom_cidx)
{
    char *BASE_ATOMS[] = { " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 ",
        " N6 ", " O2 ", " O4 ", " C5M", " O6 ", " N2 ", " N4 "
    };
    char *PRBB_ATOMS[] = { " N  ", " C  ", " O  " };
    long num_dna = sizeof BASE_ATOMS / sizeof BASE_ATOMS[0] - 1;
    long num_pro = sizeof PRBB_ATOMS / sizeof PRBB_ATOMS[0] - 1;
    long i, ib, ie, j;

    for (i = 1; i <= num_residue; i++) {
        if (res_type[i] == -2)
            continue;
        ib = seidx[i][1];
        ie = seidx[i][2];
        for (j = ib; j <= ie; j++) {
            if (atom_cidx[1][j] == 3)  /* skip H atoms */
                continue;
            if (res_type[i] == -1) {  /* amino-acid */
                if (num_strmatch(AtomName[j], PRBB_ATOMS, 0, num_pro))
                    atom_cidx[2][j] = -WITH_BKBN;
                else  /* side chain */
                    atom_cidx[2][j] = -WITH_BASE;
            } else {  /* nucleotide */
                if (num_strmatch(AtomName[j], BASE_ATOMS, 0, num_dna))
                    atom_cidx[2][j] = WITH_BASE;
                else
                    atom_cidx[2][j] = WITH_BKBN;
            }
        }
    }
}

long number_of_aa(long num_residue, long *res_type)
{
    long i, num = 0;

    for (i = 1; i <= num_residue; i++)
        if (res_type[i] == -1)
            num++;

    return num;
}

long number_of_nt(long num_residue, long *res_type)
{
    long i, num = 0;

    for (i = 1; i <= num_residue; i++)
        if (res_type[i] == 0 || res_type[i] == 1)
            num++;

    return num;
}

void get_snap_par(double **rot1, double *org1, double **rot2, double *org2,
                  char *direction, double *trs_dist, double *rot_dist, double *pars,
                  double *orgP, double **rotP)
{
    long i;
    double temp[4], xorg[4];
    double dsum = 0.0, **xmst;

    xmst = dmatrix(1, 3, 1, 3);

    strcpy(direction, "plus");
    for (i = 1; i <= 3; i++)  /* dot(z1, z2) */
        dsum += rot1[i][3] * rot2[i][3];

    if (dsum < 0.0) {  /* reverse y- and z-axes of frame #2 */
        reverse_y_z_columns(rot2);
        strcpy(direction, "minus");
    }

    bpstep_par(rot1, org1, rot2, org2, pars, xmst, xorg);
    sgl_helix(rot1, rot2, rot_dist, temp);  /* rotation axis 'temp' not used */

    ddxyz(org1, org2, temp);  /* frame 1 as reference */
    multi_vec_matrix(temp, 3, rot1, 3, 3, orgP);

    *trs_dist = veclen(temp);  /* same as veclen(orgP) */
    if (orgP[3] < 0.0)  /* aa is BELOW the mean bp plane */
        *trs_dist = -*trs_dist;

    transpose_matrix(rot1, 3, 3, xmst);
    multi_matrix(xmst, 3, 3, rot2, 3, 3, rotP);

    free_dmatrix(xmst, 1, 3, 1, 3);
}

void write_snap_par(FILE * fp, char *direction, double dist, double rot_ang,
                    double *pars, double *orgP, double **rotP)
{
    double ab_angle, x[4], o[4];
    long i, j;

    for (i = 1; i <= 3; i++) {
        x[i] = rotP[i][1];  /* C-alpha ---> C-beta as the x-axis */
        o[i] = -orgP[i];  /* C-alpha ---> bp origin */
    }
    ab_angle = magang(x, o);

    fprintf(fp, "\t%s %10.3f %10.3f %10.3f\n", direction, dist, rot_ang, ab_angle);
    for (j = 1; j <= 6; j++)
        fprintf(fp, "%10.3f", pars[j]);
    fprintf(fp, "\n");

    for (j = 1; j <= 3; j++)
        fprintf(fp, "%10.3f", orgP[j]);
    fprintf(fp, "\n");

    for (i = 1; i <= 3; i++) {
        for (j = 1; j <= 3; j++)
            fprintf(fp, "%10.3f", rotP[j][i]);
        fprintf(fp, "\n");
    }
}

void write_atom_xyz(FILE * fp, char *fmt, double *xyz, double dft)
{
    long i;

    for (i = 1; i <= 3; i++)
        fprintf(fp, fmt, (xyz == NULL) ? dft : xyz[i]);
    fprintf(fp, "\n");
}

long set2frame(long inum, long *ivec, long **seidx, double **xyz, double *morg,
               double **mst, long *serial, double **xyz_pair)
{
    long i, j, k, num = 0;

    for (i = 1; i <= inum; i++) {
        k = ivec[i];
        for (j = seidx[k][1]; j <= seidx[k][2]; j++) {
            num++;
            cpxyz(xyz[j], xyz_pair[num]);
            serial[num] = j;  /* index for extracting AtomName, ResName etc */
        }
    }
    change_xyz(0, morg, mst, num, xyz_pair);

    return num;
}

long set2frameCa(long inum, long *ivec, long **seidx, long *res_type,
                 char **AtomName, double **xyz, double *morg, double **mst,
                 long *serial, double **xyz_pair)
{
    long i, j, k, num = 0;

    for (i = 1; i <= inum; i++) {
        k = ivec[i];
        if (res_type[k] != -1) {
            for (j = seidx[k][1]; j <= seidx[k][2]; j++) {
                num++;
                cpxyz(xyz[j], xyz_pair[num]);
                serial[num] = j;  /* index for extracting AtomName, ResName etc */
            }
        } else {  /* amino-acid, select only Ca atom */
            for (j = seidx[k][1]; j <= seidx[k][2]; j++)
                if (strcmp(AtomName[j], " CA ") == 0) {
                    num++;
                    cpxyz(xyz[j], xyz_pair[num]);
                    serial[num] = j;
                    break;
                }
        }
    }
    change_xyz(0, morg, mst, num, xyz_pair);

    return num;
}

long get_nextModelNumber(char *pdbfile)
{
    char *p0;
    long k, num = 0;
    FILE *fp;

    if (exist_file(pdbfile)) {
        fp = open_file(pdbfile, "r");
        while ((p0 = my_getline(fp)) != NULL) {
            if (str_pmatch(p0, "MODEL ") && (sscanf(p0, "%*s %ld", &k) == 1))
                num++;
            if (num != k)
                fprintf(stderr, "model #s not consecutive [%ld vs %ld]\n", num, k);
            free(p0);
        }
        close_file(fp);
    }

    return num + 1;  /* the new model number */
}

void output_naa_str(char *filename, char *fmode, char *idmsg, long num, long *serial,
                    char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                    double **xyz_pair, char **Miscs, long out_org)
/* output the superimposed structures of interacting nt and aa */
{
    char *p;
    long i, k, num_str;
    FILE *fp;

    fp = open_file(filename, fmode);

    if (!is_empty_string(idmsg)) {
        num_str = get_nextModelNumber(filename);
        fprintf(fp, "%6s    %4ld\n", "MODEL ", num_str);
        fprintf(fp, "REMARK    Section #%4.4ld\n", num_str);
        fprintf(fp, "REMARK    %s\n", idmsg);
    }

    for (i = 1; i <= num; i++) {
        k = serial[i];
        p = Miscs[k];
        fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n", (p[0] == 'A')
                ? "ATOM  " : "HETATM", i, AtomName[k], p[1], ResName[k], ChainID[k],
                ResSeq[k], p[2], xyz_pair[i][1], xyz_pair[i][2], xyz_pair[i][3], p + 3);
    }

    if (out_org)
        fprintf(fp, "HETATM 9999 XX   ORG X 999       0.000   0.000   0.000\n");

    if (!is_empty_string(idmsg))
        fprintf(fp, "ENDMDL\n");
    else
        fprintf(fp, "END\n");

    close_file(fp);
}

void cleanup_files(long renew, long cleanup)
{
    static char *SAA[] = { AA_LIST };
    static char *RY[] = { "AT", "GC" };
    static char *DNA = DNA_BASE;
    static char *ALC[] = { "atom_lkg", "ablk_lkg", "bblk_lkg", "pblk_lkg",
        "pall_lkg", "snap_all", "snap_blk"
    };
    static char *MISC[] = { "snap_bp.dat", "snap_nts.pdb", "snap_options" };
    char filename[BUF512], suffix[BUF512];
    long i, j, k, num_saa, num_ry, num_dna, num_alc, num_misc;

    if (!renew && !cleanup)
        return;

    num_saa = sizeof SAA / sizeof SAA[0] - 1;
    num_ry = sizeof RY / sizeof RY[0] - 1;
    num_dna = strlen(DNA) - 1;
    num_alc = sizeof ALC / sizeof ALC[0] - 1;
    num_misc = sizeof MISC / sizeof MISC[0] - 1;

    for (i = 0; i <= num_dna; i++) {
        for (j = 0; j <= num_saa; j++) {
            for (k = -1; k <= 3; k++) {
                if (k < 0)  /* for backward compatibility */
                    strcpy(suffix, "");
                else
                    sprintf(suffix, "_%ld", k);
                sprintf(filename, "%c-%s%s.par", DNA[i], SAA[j], suffix);
                remove_file(filename);
                sprintf(filename, "%c-%s%s.pdb", DNA[i], SAA[j], suffix);
                remove_file(filename);
            }
        }
    }

    for (i = 0; i <= num_ry; i++) {
        for (j = 0; j <= num_saa; j++) {
            for (k = -1; k <= 3; k++) {
                if (k < 0)  /* for backward compatibility */
                    strcpy(suffix, "");
                else
                    sprintf(suffix, "_%ld", k);
                sprintf(filename, "%s-%s%s.par", RY[i], SAA[j], suffix);
                remove_file(filename);
                sprintf(filename, "%s-%s%s.pdb", RY[i], SAA[j], suffix);
                remove_file(filename);
            }
        }
    }

    for (i = 0; i <= num_alc; i++) {
        sprintf(filename, "%s.alc", ALC[i]);
        remove_file(filename);
    }

    for (i = 0; i <= num_misc; i++)
        remove_file(MISC[i]);

    if (cleanup)
        fatal("done with cleaning up snap/poco files.\n");
}
