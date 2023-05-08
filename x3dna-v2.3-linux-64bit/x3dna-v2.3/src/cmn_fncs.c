#include "x3dna.h"

long set_3letter_base_pdb(char *res_name, char *spdb)
{
    char str[BUF32];

    strcpy(str, res_name);
    upperstr(str);
    cvtstr_c1toc2(str, ' ', '_');

    sprintf(spdb, "Atomic_%s.pdb", str);
    if (exist_file(spdb))
        return TRUE;

    sprintf(spdb, "%sconfig/Atomic_%s.pdb", Gvars.X3DNA_HOMEDIR, str);

    return exist_file(spdb);
}

void set_std_base_pdb(char *bdir, long irna, char bname, char *spdb)
{
    char mb[BUF32], str[BUF512];

    UNUSED_PARAMETER(bdir);  /* for compatibility */

    strcpy(mb, irna ? "r" : "");

    /* current directory first */
    if (isupper((int) bname))
        sprintf(spdb, "%sAtomic_%c.pdb", mb, bname);
    else  /* for modified bases */
        sprintf(spdb, "%sAtomic.%c.pdb", mb, bname);
    if (exist_file(spdb))
        return;

    /* the system directory */
    sprintf(str, "%sconfig/%s", Gvars.X3DNA_HOMEDIR, spdb);
    strcpy(spdb, str);
}

void set_std_base_pdb00(char *bdir, long irna, char bname, char *spdb)
{
    if (isupper((int) bname))
        sprintf(spdb, "%s%sAtomic_%c.pdb", bdir, irna ? "r" : "", bname);
    else  /* for modified bases */
        sprintf(spdb, "%s%sAtomic.%c.pdb", bdir, irna ? "r" : "", bname);
}

void print_used_time(time_t time0)
{
    char str[BUF512];
    double dtime;
    long minute_secs = 60, hour_secs = 60 * 60, day_secs = 24 * 60 * 60;
    long days, hours, minutes, seconds;

    dtime = difftime(time(NULL), time0);
    sprintf(str, "%.0f", dtime);
    if (sscanf(str, "%ld", &seconds) != 1)
        fatal("wrong time format\n");

    days = seconds / day_secs;
    seconds %= day_secs;

    hours = seconds / hour_secs;
    seconds %= hour_secs;

    minutes = seconds / minute_secs;
    seconds %= minute_secs;

    fprintf(stderr, "\nTime used: %2.2ld:%2.2ld:%2.2ld:%2.2ld\n", days, hours, minutes, seconds);
}

void parcat(char *str, double par, char *format, char *bstr)
{
    char temp[BUF512];

    if (par > EMPTY_CRITERION) {
        sprintf(temp, format, par);
        strcat(str, temp);
    } else
        strcat(str, bstr);
}

void print_bp_crit(miscPars * misc_pars, FILE * fp)
{
    fprintf(fp, "Base-pair criteria used: ");
    fprintf(fp, "%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f [%s]\n",
            misc_pars->hb_dist1, misc_pars->hb_dist2, misc_pars->max_dorg,
            misc_pars->max_dv, misc_pars->max_plane_angle, misc_pars->min_dNN,
            misc_pars->helix_break, misc_pars->hb_atoms);
}

/* endofline: check for and consume \r, \n, \r\n, or EOF */
static int endofline(FILE * fp, int c)
{
    int eol;

    eol = (c == '\r' || c == '\n');
    if (c == '\r') {
        c = getc(fp);
        if (c != '\n' && c != EOF)
            ungetc(c, fp);  /* read too far; put c back */
    }
    return eol;
}

static char *enlarge_cline(long *maxline, char *line)
/* for a 0-index char-array */
{
    char *newline;

    *maxline *= SFACTOR;

    if ((newline = (char *) realloc(line, (*maxline) * sizeof(char))) == NULL)
        fatal("realloc failure in enlarge_cline()\n");

    return newline;
}

char *my_getline(FILE * fp)
/* read line of arbitrary length from fp */
{
    int c;
    long i, maxline = BUF512;
    char *line = NULL;

    line = (char *) malloc(maxline * sizeof(char));
    if (line == NULL)
        fatal("out of memory in my_getline()\n");

    for (i = 0; (c = getc(fp)) != EOF && !endofline(fp, c); i++) {
        if (i >= maxline - 1)
            line = enlarge_cline(&maxline, line);
        line[i] = c;
    }
    line[i] = '\0';

    if (c == EOF && i == 0) {
        free(line);  /* to avoid memory leak here */
        line = NULL;
    }

    return line;
}

/* strtok is GREEDY, meaning empty field will be ignored:
 * e.g. \t\t. This function is created to fix the problem.
 * Note that item[] is 1-indexed*/
long csplit(char *str, char *item[], long itemsize, char sepc)
{
    char *p0, *p;
    long nitem = 0;

    if (str[0] == '\0')
        return nitem;

    p = str;
    while (*p) {  /* change '"' to ' ' */
        if (*p == '"')
            *p = ' ';
        p++;
    }

    p0 = str;
    while ((p = strchr(p0, sepc)) != NULL) {
        *p = '\0';
        item[++nitem] = trim(p0);
        if (nitem >= itemsize)
            return itemsize;
        p0 = p + 1;
    }
    item[++nitem] = trim(p0);
    return nitem;
}

/* trim leading and trailing white spaces */
char *trim(char *a)
{
    int c;

    while (isspace(c = *a))
        a++;
    for (c = strlen(a) - 1; c >= 0; c--)
        if (!isspace((int) a[c]))
            break;
    a[c + 1] = '\0';

    return a;
}

/* trim leading white spaces */
char *ltrim(char *a)
{
    int c;

    while (isspace(c = *a))
        a++;

    return a;
}

/* trim trailing white spaces */
char *rtrim(char *a)
{
    int c;

    for (c = strlen(a) - 1; c >= 0; c--)
        if (!isspace((int) a[c]))
            break;
    a[c + 1] = '\0';

    return a;
}

/* itemize a string into its components */
long itemize(char *str, char *item[], long itemsize)
{
    static char *sep_chars = " \t\r\n";
    char *p;
    long nitem = 0;

    for (p = strtok(str, sep_chars); p != NULL; p = strtok(NULL, sep_chars))
        if (nitem > itemsize) {
            fprintf(stderr, "ignoring items following %ldth\n", nitem);
            return itemsize;
        } else
            item[nitem++] = p;
    if (!nitem)
        fatal("empty input line not allowed\n");
    return nitem - 1;  /* base-pair does not count */
}

/* itemize a string into tokens, and return number of token fields: 1-index */
long item_list(char *str, char *item[], long itemsize, char *sep_chars)
{
    char *p;
    long nitem = 0;

    for (p = strtok(str, sep_chars); p != NULL; p = strtok(NULL, sep_chars)) {
        item[++nitem] = trim(p);  /* get rid of leading and trailing spaces */
        if (nitem >= itemsize)
            return itemsize;
    }

    return nitem;
}

/* right-left base in a pair */
void refs_right_left(long bnum, double **orien, double **org, double **r1, double *o1,
                     double **r2, double *o2)
{
    long ioffset3, ioffset9;

    ioffset3 = (bnum - 1) * 3;
    ioffset9 = (bnum - 1) * 9;

    cpxyz(org[2] + ioffset3, o1);
    cpxyz(org[1] + ioffset3, o2);
    orien2mst(orien[2], ioffset9, r1);
    orien2mst(orien[1], ioffset9, r2);
}

/* base (-pair) step */
void refs_i_j(long b1, long b2, double *bp_orien, double *bp_org, double **r1,
              double *o1, double **r2, double *o2)
{
    ref_frame_i(b1, bp_orien, bp_org, r1, o1);
    ref_frame_i(b2, bp_orien, bp_org, r2, o2);
}

/* get the reference frame for "bnum" */
void ref_frame_i(long bnum, double *bp_orien, double *bp_org, double **r, double *o)
{
    long ioffset3, ioffset9;

    ioffset3 = (bnum - 1) * 3;
    ioffset9 = (bnum - 1) * 9;

    cpxyz(bp_org + ioffset3, o);
    orien2mst(bp_orien, ioffset9, r);
}

/* set orientation vector with ioffset from a 3-by-3 column-wise matrix */
void mst2orien(double *orien_vec, long ioffset, double **mst)
{
    long i, ik, j;

    for (i = 1; i <= 3; i++) {
        ik = ioffset + (i - 1) * 3;
        for (j = 1; j <= 3; j++)
            orien_vec[ik + j] = mst[j][i];
    }
}

/* get 3-by-3 column-wise matrix from orientation vector with ioffset */
void orien2mst(double *orien_vec, long ioffset, double **mst)
{
    long i, ik, j;

    for (i = 1; i <= 3; i++) {
        ik = ioffset + (i - 1) * 3;
        for (j = 1; j <= 3; j++)
            mst[j][i] = orien_vec[ik + j];
    }
}

/* given 3 vectors (x, y & z), get matrix (column-wise) */
void x_y_z_2_mtx(double *x, double *y, double *z, double **mtx)
{
    long i;

    for (i = 1; i <= 3; i++) {
        mtx[i][1] = x[i];
        mtx[i][2] = y[i];
        mtx[i][3] = z[i];
    }
}

/* given a matrix (column-wise), get 3 vectors (x, y & z) */
void mtx_2_x_y_z(double **mtx, double *x, double *y, double *z)
{
    long i;

    for (i = 1; i <= 3; i++) {
        x[i] = mtx[i][1];
        y[i] = mtx[i][2];
        z[i] = mtx[i][3];
    }
}

/* Averaging method 1 (CEHS way): always used for a single base-pair */
void cehs_average(long inum_base, long *ivec, double **orien, double **org, double **mst,
                  double *morg)
{
    double pars[7], orgn[4], **bi, **mstn;
    long i, ik;

    mstn = dmatrix(1, 3, 1, 3);
    bi = dmatrix(1, 3, 1, 3);

    ik = ivec[1];  /* middle frame initialized to 1st base */
    cpxyz(org[ik], morg);
    orien2mst(orien[ik], 0, mst);
    for (i = 1; i <= inum_base; i++) {  /* for each residue */
        ik = ivec[i];
        orien2mst(orien[ik], 0, bi);
        if (dot(&orien[ik][6], &orien[ivec[1]][6]) < 0.0)
            reverse_y_z_columns(bi);
        bpstep_par(bi, org[ik], mst, morg, pars, mstn, orgn);
        copy_dmatrix(mstn, 3, 3, mst);
        cpxyz(orgn, morg);
    }

    free_dmatrix(mstn, 1, 3, 1, 3);
    free_dmatrix(bi, 1, 3, 1, 3);
}

/* Geometrical averaging (z-, x-, then y-axes): always used for > 2 bases */
void geom_average(long inum_base, long *ivec, double **orien, double **org, double **mst,
                  double *morg)
{
    long ap, i, ik, j;
    double sx[4], sy[4], sz[4];

    for (i = 1; i <= 3; i++)
        morg[i] = sx[i] = sz[i] = 0.0;  /* y-axis follows x- and z-axes */

    for (i = 1; i <= inum_base; i++) {  /* for each residue */
        ik = ivec[i];
        ap = dot(&orien[ik][6], &orien[ivec[1]][6]) < 0.0;
        for (j = 1; j <= 3; j++) {
            morg[j] += org[ik][j];
            sx[j] += orien[ik][j];  /* x-axis */
            sz[j] += (ap) ? -orien[ik][6 + j] : orien[ik][6 + j];
        }
    }
    vec_norm(sz);
    vec_orth(sx, sz);  /* x-axis normalized */
    cross(sz, sx, sy);

    x_y_z_2_mtx(sx, sy, sz, mst);
    for (i = 1; i <= 3; i++)
        morg[i] /= inum_base;
}

/* Write a base-pair or a multiplet with reference to its middle-frame:
 * There are 2 methods to calculate the middle-frame with > 2 bases:
 * [1] Get the middle-frame between 1 & 2 (m1_2) in the usual way
 *     then the middle-frame between m1_2 & 3 (m1_2_3), etc.
 *     The final result DEPENDS ON the order, i.e., m1_3_2 is different
 *     from m1_2_3, both in position and in orientation!
 * [2] Get the mean z-axis of ALL BASES as the middle-frame z-axis, and
 *     similarly for the x-axis but with a correction for orthogonality.
 *     y-axis follows a right-handed rule. Applying y-first will give
 *     slightly different result; for pan-anti tetrad which is symmetrical,
 *     averaging will give near zero values: ===> big error!
 * So [1] is used all the time ...... */
void pair2mst(long inum_base, long *ivec, char **AtomName, char **ResName, char *ChainID,
              long *ResSeq, char **Miscs, double **xyz, double **orien, double **org,
              long **seidx, double *mst_orien, double *mst_org, long **htm_water,
              miscPars * misc_pars, FILE * fp)
{
    double morg[4], **mst, **xyz_residue;
    long i, ik, j, m, tnum_res, inum = 0;
    long ivec2[BUF512];

    tnum_res = attached_residues(inum_base, ivec, ivec2, seidx, xyz, htm_water, misc_pars);

    mst = dmatrix(1, 3, 1, 3);
    init_dvector(morg, 1, 3, 0.0);
    cehs_average(inum_base, ivec, orien, org, mst, morg);

    if (mst_orien != NULL) {  /* return middle frame */
        cpxyz(morg, mst_org);
        mst2orien(mst_orien, 0, mst);
    }

    xyz_residue = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
    for (i = 1; i <= tnum_res; i++) {  /* for each residue + attachments */
        ik = ivec2[i];
        for (j = seidx[ik][1]; j <= seidx[ik][2]; j++) {
            m = j - seidx[ik][1] + 1;
            cpxyz(xyz[j], xyz_residue[m]);
        }
        if (!Gvars.ORIGINAL_COORDINATE)
            change_xyz(0, morg, mst, seidx[ik][2] - seidx[ik][1] + 1, xyz_residue);
        pdb_record(seidx[ik][1], seidx[ik][2], &inum, 1, AtomName, ResName,
                   ChainID, ResSeq, xyz_residue, Miscs, fp);
    }

    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(xyz_residue, 1, NUM_RESIDUE_ATOMS, 1, 3);
}

/* get the chi torsion angle for each base residue */
void get_chi_angle(long num_residue, long *RY, char *bseq, long **seidx, double **xyz,
                   char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                   char **Miscs, double *chi, long **idxCN)
{
    char idmsg[BUF512];
    long i, ib, ie, idx[5], j;
    double **xyz4;

    xyz4 = dmatrix(1, 4, 1, 3);
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] < 0)
            continue;  /* no base residue */
        ib = seidx[i][1];
        ie = seidx[i][2];
        get_idmsg(ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2], idmsg);
        idx[1] = find_1st_atom(" O4'", AtomName, ib, ie, idmsg);
        idx[2] = find_1st_atom(" C1'", AtomName, ib, ie, idmsg);
        if (RY[i] == 1) {
            idx[3] = find_1st_atom(" N9 ", AtomName, ib, ie, idmsg);
            idx[4] = find_1st_atom(" C4 ", AtomName, ib, ie, idmsg);
        } else {
            if (bseq[i] == 'P' || bseq[i] == 'p') {
                idx[3] = find_1st_atom(" C5 ", AtomName, ib, ie, idmsg);
                idx[4] = find_1st_atom(" C4 ", AtomName, ib, ie, idmsg);
            } else {
                idx[3] = find_1st_atom(" N1 ", AtomName, ib, ie, idmsg);
                idx[4] = find_1st_atom(" C2 ", AtomName, ib, ie, idmsg);
            }
        }
        if (idxCN != NULL) {  /* get C1* and RN9/YN1 indices */
            idxCN[i][1] = idx[2];
            idxCN[i][2] = idx[3];
        }
        for (j = 1; j <= 4; j++) {
            if (!idx[j])
                break;
            cpxyz(xyz[idx[j]], xyz4[j]);
        }
        chi[i] = (j > 4) ? torsion(xyz4) : EMPTY_NUMBER;
    }
    free_dmatrix(xyz4, 1, 4, 1, 3);
}

FILE *open_tmpfile(void)
{
    FILE *fp;

    errno = 0;
    fp = tmpfile();
    if (fp == NULL)
        fatal("open_tmpfile() failed: %s\n", strerror(errno));

    return fp;
}

FILE *open_file(char *filename, char *filemode)
{
    FILE *fp;

    errno = 0;
    if (filename == NULL)
        filename = "\0";
    if (!strcmp(filename, "stdin"))
        fp = stdin;
    else if (!strcmp(filename, "stdout"))
        fp = stdout;
    else if (!strcmp(filename, "stderr"))
        fp = stderr;
    else {
        fp = fopen(filename, filemode);
        if (fp == NULL)
            fatal("open_file <%s> failed: %s\n", filename, strerror(errno));
    }
    return fp;
}

long close_file(FILE * fp)
{
    long i;

    if (fp == NULL || fp == stdin || fp == stdout || fp == stderr)
        return 0;
    errno = 0;
    i = fclose(fp);
    if (i == EOF)
        fatal("close_file failed: %s\n", strerror(errno));
    return i;
}

long exist_file(char *filename)
{
    long status;
    FILE *fp;

    fp = fopen(filename, "r");
    status = (fp != NULL) ? 1 : 0;
    close_file(fp);
    return status;
}

void remove_file(char *filename)
{
    if (!exist_file(filename))
        return;

    if (remove(filename))
        fatal("can not remove file: %s\n", filename);
}

void rename_file(char *src, char *dst)
{
    if (!exist_file(src))
        fatal("file to be renamed <%s> does NOT exist\n", src);

    if (rename(src, dst))
        fatal("can not rename file <%s> to <%s>\n", src, dst);
}

void copy_file_pointer(FILE * fpi, FILE * fpo, char *msg)
{
    char str[BUF512];
    size_t num_bytes;

    while (!feof(fpi)) {
        num_bytes = fread(str, 1, BUF512, fpi);
        if (num_bytes == 0)
            fprintf(stderr, "zero bytes? [%s]\n", msg);
        if (fwrite(str, 1, num_bytes, fpo) != num_bytes)
            fatal("file %s error\n", msg);
    }
}

/* copy file contents from src to dst */
void cpcat_file(char *src, char *dst, char *method)
{
    FILE *fpi, *fpo;

    if (!strcmp(src, dst)) {
        fprintf(stderr, "same source/destination file name: <%s>\n", src);
        return;
    }

    fpi = open_file(src, "rb");

    if (strstr("catenation", method) || strstr("append", method)) {
        fpo = open_file(dst, "ab");
        copy_file_pointer(fpi, fpo, "concatenation");
        close_file(fpo);

    } else {
        fpo = open_file(dst, "wb");
        copy_file_pointer(fpi, fpo, "copy");
        close_file(fpo);
    }

    close_file(fpi);
}

/* change to upper case, and return string length */
long upperstr(char *a)
{
    long nlen = 0;

    while (*a) {
        nlen++;
        if (islower((int) *a))
            *a = toupper((int) *a);
        a++;
    }
    return nlen;
}

/* change to lower case, and return string length */
long lowerstr(char *a)
{
    long nlen = 0;

    while (*a) {
        nlen++;
        if (isupper((int) *a))
            *a = tolower((int) *a);
        a++;
    }
    return nlen;
}

char *my_strdup(const char *src)
{
    char *dst;

    dst = cvector(0, strlen(src));
    strcpy(dst, src);

    return dst;
}

/* print char 'x' n-times to stream fp */
void print_sep(FILE * fp, char x, long n)
{
    long i;

    for (i = 1; i <= n; i++)
        if (fputc(x, fp) == EOF)
            fatal("error writing characters to the stream\n");
    if (fputc('\n', fp) == EOF)
        fatal("error writing '\n' to the stream\n");
}

/* check if '/' exist at the end of the original BDIR */
void check_slash(char *BDIR)
{
    char *pchar;
    long n;

    pchar = strrchr(BDIR, '/');
    n = strlen(BDIR);
    if (pchar && (pchar - BDIR != n - 1)) {
        BDIR[n] = '/';
        BDIR[n + 1] = '\0';
    }
}

/* delete '/' at the end of the str */
void delete_end_slash(char *str)
{
    char *pchar;
    long n;

    pchar = strrchr(str, '/');
    n = strlen(str);
    if (pchar - str == n - 1)
        str[n - 1] = '\0';
}

/* return a pointer 1-char following the last slash or to str w/o '/' */
char *basename(char *str)
{
    char *p_lslash;

    p_lslash = strrchr(str, '/');
    if (p_lslash == NULL)
        return str;
    else
        return p_lslash + 1;
}

long lround(double d)
{
    return (long) ((d > 0.0) ? d + 0.5 : d - 0.5);
}

/* get rid of the extension in a file name */
void del_extension(char *fullname, char *okname)
{
    char *pchar, *bname;
    size_t i;

    bname = basename(fullname);

    pchar = strrchr(bname, '.');
    if (pchar == NULL)
        strcpy(okname, bname);
    else {
        i = pchar - bname;
        strncpy(okname, bname, i);
        okname[i] = '\0';
    }
}

/* get the base name without last extension */
void bname_noext(char *src, char *dst)
{
    char str[BUF512];

    strcpy(str, basename(src));
    del_extension(str, dst);
}

/* exit program with an error message */
void fatal(char *fmt, ...)
{
    va_list args;

    if (strlen(fmt) > 0) {
        va_start(args, fmt);
        vfprintf(stderr, fmt, args);
        va_end(args);
    }

    exit(1);
}

/* print HEADER, TITLE etc records to parameter output file for easy reference */
void print_pdb_title(char *pdbfile, char *chain_list, FILE * fp)
{
    char str[BUF512];
    static char *titles[] = { "HEADER", "TITLE ", "AUTHOR", "COMPND", "SOURCE", "EXPDTA", "KEYWDS",
        "REVDAT", "JRNL  ", "HELIX ", "SHEET ", "TURN  "
    };
    long i, num, nlen;
    FILE *fpp;

    if (!Gvars.HEADER)
        return;

    num = sizeof titles / sizeof titles[0];

    fpp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fpp) != NULL) {
        nlen = upperstr(str);
        if (!strncmp(str, "ATOM  ", 6) || !strncmp(str, "HETATM", 6) || !strncmp(str, "END", 3))
            break;
        if (Gvars.HEADER > TRUE) {
            fprintf(fp, "%s", str);
            continue;
        }

        /* default case */
        if (nlen >= 6) {  /* at least 6 characters */
            for (i = 0; i < num - 3; i++)
                if (!strncmp(str, titles[i], 6))
                    fprintf(fp, "%s", str);
            if ((!strncmp(str, "HELIX ", 6) || !strncmp(str, "TURN  ", 6)) &&
                (strchr(chain_list, '*') || strchr(chain_list, str[19])))
                fprintf(fp, "%s", str);
            if (!strncmp(str, "SHEET ", 6) &&
                (strchr(chain_list, '*') || strchr(chain_list, str[21])))
                fprintf(fp, "%s", str);
        }
    }
    close_file(fpp);
}

static long is_end_of_structure_to_process(char *str)
{
    if (str_pmatch(str, "END")) {
        if (Gvars.ALL_MODEL)
            return str_pmatch(str, "ENDMDL") ? FALSE : TRUE;
        else  /* also matches ENDMDL */
            return TRUE;
    } else
        return FALSE;
}

static double get_occupancy(long nlen, char *str, char *pdbfile)
{
    char temp[BUF512];
    double occupancy;

    if (nlen < 60 || !Gvars.OCCUPANCY)  /* no checking for occupancy */
        return 1.0;

    strncpy(temp, str + 54, 6);  /* occupancy */
    temp[6] = '\0';
    if (sscanf(temp, "%6lf", &occupancy) != 1)
        fatal("error reading occupancy in file [%s]: '%s'\n", pdbfile, str);

    return occupancy;
}

/* The last column of z-coordinate is 54. However, some non-compliant
 * PDB could be less then 54 column per ATOM/HETATM record. So here
 * use 52 to account for such special cases */
#define Zcol  52

/* number of atom records in a PDB file */
long number_of_atoms(char *pdbfile, long hetatm, char *ALT_LIST)
{
    char str[BUF512], *pchar;
    long n = 0, nlen;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        if ((pchar = strrchr(str, '\n')) != NULL)
            str[pchar - str] = '\0';
        nlen = upperstr(str);

        if (is_end_of_structure_to_process(str))
            break;

        if (nlen >= Zcol && (!strncmp(str, "ATOM", 4)
                             || (hetatm && !strncmp(str, "HETATM", 6)))
            && (strchr(ALT_LIST, '*') || strchr(ALT_LIST, str[16]))
            && get_occupancy(nlen, str, pdbfile) > 0)
            n++;
    }
    close_file(fp);

    return n;
}

/* read in a PDB file and do some processing
 * Miscs[][NMISC]: H/A, altLoc, iCode, occ./tempFac./segID/element/charge
 *           col#   0     1       2       3-28 [combined together]
 * 29-32, 4 characters for Model Numbers (in PDB MODEL record, col: 11-14)
 * Element symbol: column# 25-26, right justified
 * 80 - 54 = 26 characters from ATOM/HETATM records 26 + 2 = 28
 *      set the 29th column to '\0' for end of the ATOM/HETATM record
 * 30 -- 33 for the 4-character model number; 34 is set to '\0'
 *   so we have NMISC = 34 */
long read_pdb(char *pdbfile, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs, long hetatm, char *ALT_LIST)
{
    char *p0, *pn;
    char rname[4], str[BUF512], str0[BUF512], temp[BUF512];
    double occupancy;
    long i, n = 0, nlen, k, occ_chk = FALSE, modelNum = 0;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while ((p0 = my_getline(fp)) != NULL) {
        strcpy(str, p0);
        free(p0);

        strcpy(str0, str);
        nlen = upperstr(str);
        if (!strncmp(str, "MODEL ", 6)  /* model record */
            &&(sscanf(str, "%*s %ld", &k) == 1))
            modelNum = k;

        if (is_end_of_structure_to_process(str))
            break;

        if (nlen >= Zcol &&  /* same as in number_of_atoms */
            (!strncmp(str, "ATOM", 4)
             || (hetatm && !strncmp(str, "HETATM", 6)))
            && (strchr(ALT_LIST, '*') || strchr(ALT_LIST, str[16]))) {

            occupancy = get_occupancy(nlen, str, pdbfile);
            if (occupancy <= 0) {  /* ignore 0-occupancy atom */
                if (!occ_chk) {
                    occ_chk = TRUE;
                    fprintf(stderr, "[i] File '%s' with atom occupancy <= 0 [%s]\n", pdbfile, str);
                }
                continue;
            }

            n++;

            if (AtomSNum != NULL) {  /* read in original atom serial number */
                strncpy(temp, str + 6, 5);
                temp[5] = '\0';  /* 5 digits */
                if (sscanf(temp, "%5ld", &AtomSNum[n]) != 1) {
                    fprintf(stderr, "Atom serial number ? %s ?\n", str);
                    AtomSNum[n] = n;  /* sequential number */
                }
            }

            strncpy(AtomName[n], str + 12, 4);
            AtomName[n][4] = '\0';

            if (Gvars.AtomName0 && Gvars.Name0)
                strcpy(Gvars.AtomName0[n], AtomName[n]);

            /* for very special case 'P   ', 'N9  ' etc. Fei XU: also "MG  " 355d */
            pn = AtomName[n];  /* as a shorthand */
            if ((pn[0] != ' ' && !isdigit((int) pn[0]))
                && (pn[1] == ' ' || isdigit((int) pn[1])) && pn[3] == ' ') {
                strncpy(temp, pn, 3);
                temp[3] = '\0';
                sprintf(pn, " %s", temp);
            } else if (pn[0] == ' ' && pn[1] == ' ' && isdigit((int) pn[3])) {  /* as in "  N1" */
                strncpy(temp, pn + 2, 2);
                temp[2] = '\0';
                sprintf(pn, " %s ", temp);
            } else if (is_equal_string(pn, "   P") || is_equal_string(pn, "P   "))
                strcpy(pn, " P  ");
            else if (is_equal_string(pn, "OP1 "))
                strcpy(pn, " O1P");
            else if (is_equal_string(pn, "OP2 "))
                strcpy(pn, " O2P");

            if (AtomName[n][3] == '*')  /* * to ' */
                AtomName[n][3] = '\'';
            if (!strcmp(AtomName[n], " O1'"))  /* O1' to O4' */
                strcpy(AtomName[n], " O4'");
            if (!strcmp(AtomName[n], " OL ") || !strcmp(AtomName[n], " OP1"))  /* OL/OP1 to O1P */
                strcpy(AtomName[n], " O1P");
            if (!strcmp(AtomName[n], " OR ") || !strcmp(AtomName[n], " OP2"))  /* OR/OP2 to O2P */
                strcpy(AtomName[n], " O2P");
            if (!strcmp(AtomName[n], " OP3"))  /* OP3 to O3P */
                strcpy(AtomName[n], " O3P");
            if (!strcmp(AtomName[n], " C5A"))  /* C5A to C5M */
                strcpy(AtomName[n], " C5M");
            if (!strcmp(AtomName[n], " O5T"))  /* terminal O5' */
                strcpy(AtomName[n], " O5'");
            if (!strcmp(AtomName[n], " O3T"))  /* terminal O3' */
                strcpy(AtomName[n], " O3'");

            strncpy(rname, str + 17, 3);  /* residue name */
            rname[3] = '\0';
            if (Gvars.ResName0 && Gvars.Name0)
                strcpy(Gvars.ResName0[n], rname);

            /* delete ending spaces as in "C  " */
            for (i = 1; i <= 2; i++)
                if (rname[2] == ' ') {
                    rname[2] = rname[1];
                    rname[1] = rname[0];
                    rname[0] = ' ';
                }
            if (rname[2] == ' ')
                fatal("==> residue name [%s] field empty <==\n", str);
            strcpy(ResName[n], rname);

            ChainID[n] = Gvars.CHAIN_CASE ? str0[21] : str[21];
            strncpy(temp, str + 22, 4);
            temp[4] = '\0';
            if (sscanf(temp, "%4ld", &ResSeq[n]) != 1) {
                fprintf(stderr, "residue #? ==> %.54s\n", str);
                ResSeq[n] = 9999;
            }
            strncpy(temp, str + 30, 8);
            temp[8] = '\0';
            if (sscanf(temp, "%8lf", &xyz[n][1]) != 1)
                fatal("error reading x-coordinate\n");
            strncpy(temp, str + 38, 8);
            temp[8] = '\0';
            if (sscanf(temp, "%8lf", &xyz[n][2]) != 1)
                fatal("error reading y-coordinate\n");
            strncpy(temp, str + 46, 8);
            temp[8] = '\0';
            if (sscanf(temp, "%8lf", &xyz[n][3]) != 1)
                fatal("error reading z-coordinate\n");
            if (Miscs != NULL) {
                Miscs[n][0] = str[0];  /* H for HETATM, A for ATOM */
                Miscs[n][1] = str[16];  /* alternative location indicator */
                Miscs[n][2] = str[26];  /* code of insertion residues */
                if (nlen >= 54)
                    strncpy(Miscs[n] + 3, str + 54, 26);  /* upto the 80th column */
                Miscs[n][29] = '\0';
                sprintf(Miscs[n] + 30, "%4.4ld", modelNum);
                Miscs[n][NMISC] = '\0';  /* just to make sure */
            }
        }
    }
    close_file(fp);

    if (!n)
        fprintf(stderr, "PDB file <%s> has NO ATOM/HETATM records\n", pdbfile);

    return n;
}

/* free all PDB relevant arrays to make code clean and short */
void free_pdb(long num, long *AtomSNum, char **AtomName, char **ResName, char *ChainID,
              long *ResSeq, double **xyz, char **Miscs)
{
    if (AtomSNum != NULL)
        free_lvector(AtomSNum, 1, num);
    if (AtomName != NULL)
        free_cmatrix(AtomName, 1, num, 0, 4);
    if (ResName != NULL)
        free_cmatrix(ResName, 1, num, 0, 3);
    if (ChainID != NULL)
        free_cvector(ChainID, 1, num);
    if (ResSeq != NULL)
        free_lvector(ResSeq, 1, num);
    if (xyz != NULL)
        free_dmatrix(xyz, 1, num, 1, 3);
    if (Miscs != NULL)
        free_cmatrix(Miscs, 1, num, 0, NMISC);
}

/* reset xyz coordinates for PDB and ALCHEMY formats */
void reset_xyz(long num, double **xyz, char *fmt)
{
    double ave_xyz[4], max_xyz[4], min_xyz[4];

    /* check the range of the coordinates */
    max_dmatrix(xyz, num, 3, max_xyz);
    min_dmatrix(xyz, num, 3, min_xyz);
    ave_dmatrix(xyz, num, 3, ave_xyz);

    if (max_dvector(max_xyz, 1, 3) > 9999.99) {
        fprintf(stderr, "xyz coordinate over %s limit. reset origin"
                " to geometrical center\n", fmt);
        move_position(xyz, num, 3, ave_xyz);
    } else if (min_dvector(min_xyz, 1, 3) < -999.99) {
        fprintf(stderr, "xyz coordinate under %s limit. reset origin"
                " to minimum xyz coordinates\n", fmt);
        move_position(xyz, num, 3, min_xyz);
    }
}

void deduce_misc(char **Miscs, char **AtomName, long i, char *str)
{
    static char asym[3];

    if (Miscs == NULL || is_equal_string(Miscs[i], "A  ")) {
        aname2asym(AtomName[i], asym, Gvars.NUM_SATOM, Gvars.ATOMLIST);
        /* ATOM/HETATM, altLoc, iCode, occupancy, temp-factor, unused 67..76, atom-symbol, and charge */
        sprintf(str, "A  %6.2f%6.2f          %2.2s  ", 1.0, 1.0, asym);
    } else
        strcpy(str, Miscs[i]);
}

static void cvt_3letter_nts(char *rname)
{
    if (is_equal_string(rname, "  A") || is_equal_string(rname, " DA"))
        strcpy(rname, "ADE");
    else if (is_equal_string(rname, "  C") || is_equal_string(rname, " DC"))
        strcpy(rname, "CYT");
    else if (is_equal_string(rname, "  G") || is_equal_string(rname, " DG"))
        strcpy(rname, "GUA");
    else if (is_equal_string(rname, "  T") || is_equal_string(rname, " DT"))
        strcpy(rname, "THY");
    else if (is_equal_string(rname, "  U"))
        strcpy(rname, "URA");
}

long is_dna_with_backbone(long ib, long ie, char **AtomName)
{
    long i, P = FALSE;

    for (i = ib; i <= ie; i++) {
        if (is_equal_string(AtomName[i], " O2'"))
            return FALSE;  /* taken as RNA */
        if (!P && is_equal_string(AtomName[i], " P  "))
            P = TRUE;
    }

    return P;
}

static void cvt_pdbv3_name(char *aname)
{
    if (is_equal_string(aname, " O1P"))
        strcpy(aname, " OP1");
    else if (is_equal_string(aname, " O2P"))
        strcpy(aname, " OP2");
    else if (is_equal_string(aname, " C5M"))
        strcpy(aname, " C7 ");
}

static void cvt_pdbv3_dna(char *rname)
{
    if (is_equal_string(rname, "  A"))
        strcpy(rname, " DA");
    else if (is_equal_string(rname, "  C"))
        strcpy(rname, " DC");
    else if (is_equal_string(rname, "  G"))
        strcpy(rname, " DG");
    else if (is_equal_string(rname, "  T"))
        strcpy(rname, " DT");
}

void normalize_resName_atomName(long is_dna, const char *rname0, const char *aname0,
                                char *rname, char *aname)
{
    strcpy(aname, aname0);
    strcpy(rname, rname0);

    if (Gvars.PDBV3) {
        cvt_pdbv3_name(aname);
        if (is_dna)
            cvt_pdbv3_dna(rname);

    } else if (Gvars.THREE_LETTER_NTS)
        cvt_3letter_nts(rname);
}

/* write out ATOM and HETATM record: xyz could be 1 to [ie - ib + 1] */
void pdb_record(long ib, long ie, long *inum, long idx, char **AtomName, char **ResName,
                char *ChainID, long *ResSeq, double **xyz, char **Miscs, FILE * fp)
{
    char rname[4], aname[5], str[BUF512];
    long i, j, is_dna;

    is_dna = Gvars.PDBV3 && is_dna_with_backbone(ib, ie, AtomName);

    for (i = ib; i <= ie; i++) {
        deduce_misc(Miscs, AtomName, i, str);
        if (Gvars.AtomName0) {
            strcpy(aname, Gvars.AtomName0[i]);
            strcpy(rname, Gvars.ResName0[i]);
        } else
            normalize_resName_atomName(is_dna, ResName[i], AtomName[i], rname, aname);

        j = (idx) ? i - ib + 1 : i;
        fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n", (str[0] == 'A')
                ? "ATOM  " : "HETATM", ++*inum, aname, str[1], rname, ChainID[i],
                ResSeq[i], str[2], xyz[j][1], xyz[j][2], xyz[j][3], str + 3);
    }
}

void write_pdb(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               double **xyz, char **Miscs, char *pdbfile)
{
    long inum = 0;
    FILE *fp;

    reset_xyz(num, xyz, "f8.3");

    fp = open_file(pdbfile, "w");
    fprintf(fp, "REMARK    %s\n", Gvars.X3DNA_VER);
    pdb_record(1, num, &inum, 0, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, fp);
    fprintf(fp, "END\n");
    close_file(fp);
}

/* output structure in PDBML to hold large structures */
void write_pdbml(long xml, long num, char **AtomName, char **ResName, char *ChainID,
                 long *ResSeq, double **xyz, char *pdbfile)
{
    char *p, my_asym[3];
    long i;
    FILE *fp;

    fp = open_file(pdbfile, "w");
    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n");
    fprintf(fp, "<PDBx:datablock datablockName=\"3DNA\"\n");
    fprintf(fp, "   xmlns:PDBx=\"http://deposit.pdb.org/pdbML/pdbx-v0.905.xsd\"\n");
    fprintf(fp, "   xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
    fprintf(fp, "   xsi:schemaLocation=\"http://deposit.pdb.org/pdbML/pdbx-v0.905.xsd"
            " pdbx-v0.905.xsd\">\n");
    fprintf(fp, "   <PDBx:citationCategory>\n");
    fprintf(fp, "      <PDBx:citation id=\"primary\">\n");
    fprintf(fp, "         <PDBx:title>3DNA: a software package for the analysis,"
            " rebuilding and visualization of three-dimensional nucleic acid"
            " structures.</PDBx:title>\n");
    fprintf(fp, "         <PDBx:journal_abbrev>Nucleic Acids Res.</PDBx:journal_abbrev>\n");
    fprintf(fp, "         <PDBx:journal_volume>31</PDBx:journal_volume>\n");
    fprintf(fp, "         <PDBx:page_first>5108</PDBx:page_first>\n");
    fprintf(fp, "         <PDBx:page_last>5121</PDBx:page_last>\n");
    fprintf(fp, "         <PDBx:year>2003</PDBx:year>\n");
    fprintf(fp, "         <PDBx:journal_id_ASTM>%s</PDBx:journal_id_ASTM>\n", Gvars.X3DNA_VER);
    fprintf(fp,
            "         <PDBx:pdbx_database_id_PubMed>12930962</PDBx:pdbx_database_id_PubMed>\n");
    fprintf(fp, "      </PDBx:citation>\n");
    fprintf(fp, "   </PDBx:citationCategory>\n");
    fprintf(fp, "   <PDBx:citation_authorCategory>\n");
    fprintf(fp, "      <PDBx:citation_author citation_id=\"primary\""
            " name=\"Lu, X.J.\"></PDBx:citation_author>\n");
    fprintf(fp, "      <PDBx:citation_author citation_id=\"primary\""
            " name=\"Olson, W.K.\"></PDBx:citation_author>\n");
    fprintf(fp, "   </PDBx:citation_authorCategory>\n");

    if (xml == 1)
        fprintf(fp, "   <PDBx:atom_siteCategory>\n");
    else
        fprintf(fp, "   <category_atom_record>\n");

    for (i = 1; i <= num; i++) {
        aname2asym(AtomName[i], my_asym, Gvars.NUM_SATOM, Gvars.ATOMLIST);
        p = trim(my_asym);
        if (xml == 1) {
            fprintf(fp, "      <PDBx:atom_site id=\"%ld\">\n", i);
            fprintf(fp, "         <PDBx:group_PDB>ATOM</PDBx:group_PDB>\n");
            fprintf(fp, "         <PDBx:type_symbol>%s</PDBx:type_symbol>\n", p);
            fprintf(fp, "         <PDBx:label_atom_id>%s</PDBx:label_atom_id>\n",
                    trim(AtomName[i]));
            fprintf(fp, "         <PDBx:label_comp_id>%s</PDBx:label_comp_id>\n",
                    trim(ResName[i]));
            fprintf(fp, "         <PDBx:label_asym_id>%c</PDBx:label_asym_id>\n", ChainID[i]);
            fprintf(fp, "         <PDBx:label_entity_id>1</PDBx:label_entity_id>\n");
            fprintf(fp, "         <PDBx:label_seq_id>%ld</PDBx:label_seq_id>\n", ResSeq[i]);
            fprintf(fp, "         <PDBx:Cartn_x>%.3f</PDBx:Cartn_x>\n", xyz[i][1]);
            fprintf(fp, "         <PDBx:Cartn_y>%.3f</PDBx:Cartn_y>\n", xyz[i][2]);
            fprintf(fp, "         <PDBx:Cartn_z>%.3f</PDBx:Cartn_z>\n", xyz[i][3]);
            fprintf(fp, "         <PDBx:occupancy>1.0</PDBx:occupancy>\n");
            fprintf(fp, "         <PDBx:B_iso_or_equiv>0.0</PDBx:B_iso_or_equiv>\n");
            fprintf(fp, "         <PDBx:auth_seq_id>%ld</PDBx:auth_seq_id>\n", ResSeq[i]);
            fprintf(fp, "         <PDBx:auth_comp_id>%s</PDBx:auth_comp_id>\n", trim(ResName[i]));
            fprintf(fp, "         <PDBx:auth_asym_id>%c</PDBx:auth_asym_id>\n", ChainID[i]);
            fprintf(fp, "         <PDBx:auth_atom_id>%s</PDBx:auth_atom_id>\n", trim(AtomName[i]));
            fprintf(fp, "         <PDBx:pdbx_PDB_model_num>1</PDBx:pdbx_PDB_model_num>\n");
            fprintf(fp, "      </PDBx:atom_site>\n");

        } else {
            fprintf(fp, "      <atom_record id=\"%ld\">ATOM", i);
            fprintf(fp, " 1 %c %c", ChainID[i], ChainID[i]);
            fprintf(fp, " %ld %ld", ResSeq[i], ResSeq[i]);
            fprintf(fp, " ? .");
            fprintf(fp, " %s %s", trim(ResName[i]), trim(ResName[i]));
            fprintf(fp, " %s", p);
            fprintf(fp, " %s %s", trim(AtomName[i]), trim(AtomName[i]));
            fprintf(fp, " %.3f %.3f %.3f 1.0 0.0", xyz[i][1], xyz[i][2], xyz[i][3]);
            fprintf(fp, "</atom_record>\n");
        }
    }

    if (xml == 1)
        fprintf(fp, "   </PDBx:atom_siteCategory>\n");
    else
        fprintf(fp, "   </category_atom_record>\n");
    fprintf(fp, "</PDBx:datablock>\n");
    close_file(fp);
}

void write_pdbcnt(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                  double **xyz, long **connect, char *pdbfile)
{
    long i, inum = 0, j;
    FILE *fp;

    reset_xyz(num, xyz, "f8.3");

    fp = open_file(pdbfile, "w");
    fprintf(fp, "REMARK    %s\n", Gvars.X3DNA_VER);
    pdb_record(1, num, &inum, 0, AtomName, ResName, ChainID, ResSeq, xyz, NULL, fp);

    for (i = 1; i <= num; i++)
        if (connect[i][7]) {
            fprintf(fp, "CONECT%5ld", i);
            for (j = 1; j <= connect[i][7]; j++)
                fprintf(fp, "%5ld", connect[i][j]);
            /* March-2-2005: guanming.wang.adv04@alum.dartmouth.org */
            for (j = 6 + 5 * (connect[i][7] + 1); j <= 70; j++)
                fprintf(fp, " ");
            fprintf(fp, "\n");
        }
    fprintf(fp, "END\n");

    close_file(fp);
}

void move_position(double **d, long nr, long nc, double *mpos)
{
    long i, j;

    for (i = 1; i <= nr; i++)
        for (j = 1; j <= nc; j++)
            d[i][j] -= mpos[j];
}

/* number of residues, and starting-ending indexes for each */
long **residue_idx(long num, long *ResSeq, char **Miscs, char *ChainID, char **ResName,
                   long *num_residue)
{
    char iCode;
    char **bidx;
    long i, n, **seidx, *temp;

    bidx = cmatrix(1, num, 0, 12);  /* normally 9 */
    temp = lvector(1, num);

    for (i = 1; i <= num; i++) {
        iCode = (Miscs == NULL) ? ' ' : Miscs[i][2];
        sprintf(bidx[i], "%3s%c%4ld%c", ResName[i], ChainID[i], ResSeq[i], iCode);
    }
    for (i = 1; i < num; i++)
        temp[i] = strcmp(bidx[i + 1], bidx[i]) ? 1 : 0;
    temp[num] = 1;

    n = 0;  /* get number of residues */
    for (i = 1; i <= num; i++)
        if (temp[i])
            ++n;

    seidx = lmatrix(1, n, 1, 2);  /* allocate spaces */
    n = 0;
    for (i = 1; i <= num; i++)
        if (temp[i])
            seidx[++n][2] = i;
    for (i = 2; i <= n; i++)
        seidx[i][1] = seidx[i - 1][2] + 1;
    seidx[1][1] = 1;

    *num_residue = n;

    free_cmatrix(bidx, 1, num, 0, 12);
    free_lvector(temp, 1, num);

    return seidx;
}

static void set_U_C5M(char **AtomName, double **xyz, long ib, long ie)
{
    long C5, C5M, C7;

    C5 = find_1st_atom(" C5 ", AtomName, ib, ie, "");
    C5M = find_1st_atom(" C5M", AtomName, ib, ie, "");
    C7 = find_1st_atom(" C7 ", AtomName, ib, ie, "");

    if (C5 && C7 && !C5M && (p1p2_dist(xyz[C5], xyz[C7]) < 2.0))
        strcpy(AtomName[C7], " C5M");
}

long frag_contain_metal(long ib, long ie, long *is_metal)
{
    long i;

    for (i = ib; i <= ie; i++)
        if (is_metal[i])
            return TRUE;

    return FALSE;
}

void atom_metal(long num_atoms, char **AtomName, long *is_metal)
{
    static char *metals[] = {
        "LI", "BE",
        "NA", "MG", "AL",
        " K", "CA", "SC", "TI", " V", "CR", "MN", "FE", "CO",
        "NI", "CU", "ZN", "GA",
        "RB", "SR", " Y", "ZR", "NB", "MO", "TC", "RU", "RH",
        "PD", "AG", "CD", "IN", "SN",
        "CS", "BA", "HF", "TA", " W", "RE", "OS", "IR", "PT",
        "AU", "HG", "TL", "PB", "BI",
        "FR", "RA", "RF", "DB", "SG", "BH", "HS", "MT",
        "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB",
        "DY", "HO", "ER", "TM", "YB", "LU",
        "AC", "TH", "PA", " U", "NP", "PU", "AM", "CM", "BK",
        "CF", "ES", "FM", "MD", "NO", "LR"
    };

#if 0
    static char *non_metals[] = {
        " H", "HE",
        " B", " C", " N", " O", " F", "NE",
        "SI", " P", " S", "CL", "AR",
        "GE", "AS", "SE", "BR", "KR",
        "SB", "TE", " I", "XE",
        "PO", "AT", "RN"
    };
#endif

    char atom_sym[3];
    long i, num_metal;

    num_metal = sizeof metals / sizeof metals[0] - 1;  /* minus - 1 */

    for (i = 1; i <= num_atoms; i++) {
        aname2asym(AtomName[i], atom_sym, Gvars.NUM_SATOM, Gvars.ATOMLIST);
        is_metal[i] = (num_strmatch(atom_sym, metals, 0, num_metal)) ? TRUE : FALSE;
    }
}

/* Assign residue type (taking into consideration of residues with
   only P and CA atoms [pde0128]):
   Nucleic acid purines: 1; pyrimidines: 0
   Amino acids: -1
   Ligand: -3 (at least 2 atoms in the residue & non-water, or a metal atom)
   Water: -6
   Others: -2 (1 atom in the residue & non-water) */
void residue_wtype(long num_residue, long **seidx, char **ResName, char **AtomName,
                   double **xyz, char **Miscs, long *res_type, long only_ntaa)
{
    static char *WATER[] = { WATER_LIST };
    static char *SNA[] = { NT_LIST };
    static char *SAA[] = { AA_LIST };
    long i, ib, ie, id, k, num_wat, num_sna, num_saa;
    long num_atoms, *is_metal;

    num_wat = sizeof WATER / sizeof WATER[0] - 1;
    num_sna = sizeof SNA / sizeof SNA[0] - 1;
    num_saa = sizeof SAA / sizeof SAA[0] - 1;

    num_atoms = seidx[num_residue][2];
    is_metal = lvector(1, num_atoms);
    atom_metal(num_atoms, AtomName, is_metal);

    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];
        res_type[i] = residue_ident(AtomName, xyz, Miscs, ib, ie);

        if (only_ntaa || res_type[i] != -2)
            continue;

        id = ie - ib + 1;  /* # of atoms in this residue */
        if ((id < 6 &&  /* 6-membered ring */
             (num_strmatch(" P  ", AtomName, ib, ie) ||
              num_strmatch(" C1'", AtomName, ib, ie))) ||
            num_strmatch(ResName[ib], SNA, 0, num_sna)) {
            res_type[i] = 6;  /* only backbone atoms P or C1' */
            continue;
        }

        k = find_1st_atom(" CA ", AtomName, ib, ie, "");
        if ((id < 4 && k && Miscs[k][0] == 'A') ||  /* contains CA in an ATOM record */
            num_strmatch(ResName[ib], SAA, 0, num_saa)) {
            res_type[i] = -1;
            continue;
        }

        if (num_strmatch(ResName[ib], WATER, 0, num_wat) ||
            (id == 1 && (!strcmp(AtomName[ib], " O  ") || !strcmp(AtomName[ib], " OW ")))) {
            res_type[i] = -6;
            continue;
        }

        for (k = ib; k <= ie; k++)
            if (Miscs[k][0] == 'A')  /* ATOM record */
                break;
        if (k > ie && (id >= 2 || (id == 1 && is_metal[ib])))
            res_type[i] = -3;  /* ligand: HETATM at least 2 atoms or one metal atom */
    }

    free_lvector(is_metal, 1, num_atoms);
}

/* 2o8b_C30F has RMSD 0.24; normally < 0.1 */
static double check_nt_type_by_rmsd(long *idx, long C1_prime, double **xyz)
{
    /* idx[] RA_LIST: " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 " */
    static double xyz_ring[][4] = {  /* based on Atomic_G.pdb */
        {XBIG, -1.265, 3.177, 0.000},  /* C4 -- 0 */
        {XBIG, -2.342, 2.364, 0.001},  /* N3 -- 1 */
        {XBIG, -1.999, 1.087, 0.000},  /* C2 -- 2 */
        {XBIG, -0.700, 0.641, 0.000},  /* N1 -- 3 */
        {XBIG, 0.424, 1.460, 0.000},  /* C6 -- 4 */
        {XBIG, 0.071, 2.833, 0.000},  /* C5 -- 5 */
        {XBIG, 0.870, 3.969, 0.000},  /* N7 -- 6 */
        {XBIG, 0.023, 4.962, 0.000},  /* C8 -- 7 */
        {XBIG, -1.289, 4.551, 0.000}  /* N9 -- 8 */
    };
    double **xyz1, **xyz2, **fitted_xyz, **R;
    double rmsd, org[4];
    long i, k = 0, nN = 0, num = 9;

    /* large enough to hold everything: 3 * 9 + 3 = 30 */
    xyz1 = dmatrix(1, BUF32, 1, 3);
    xyz2 = &xyz1[num + 1];
    fitted_xyz = &xyz1[2 * num + 1];
    R = &xyz1[3 * num + 1];

    for (i = 0; i < num; i++) {
        if (!idx[i])
            continue;
        k++;
        if (i == 1 || i == 3 || i == 6 || i == 8)
            nN++;

        cpxyz(xyz[idx[i]], xyz1[k]);
        cpxyz(xyz_ring[i], xyz2[k]);
    }

    if (!nN && !C1_prime)
        rmsd = DUMMY;
    else
        rmsd = ls_fitting(xyz1, xyz2, k, fitted_xyz, R, org);

    free_dmatrix(xyz1, 1, DUMMY, 1, DUMMY);

    return rmsd;
}

/* identifying a residue as follows:
 *  R-base  Y-base  amino-acid, others [default]
 *   +1        0        -1        -2 [default] */
long residue_ident(char **AtomName, double **xyz, char **Miscs, long ib, long ie)
{
    static char *RingAtom[] = { RA_LIST };
    double rmsd = XBIG, dcrt = 2.0;
    long i, n, k = 0, kr = 0, num = 9, idx[16] = { 0 };
    long CA, C, C1_prime;

    for (i = 0; i < num; i++) {  /* nine atoms */
        n = find_1st_atom(RingAtom[i], AtomName, ib, ie, "");
        if (n) {
            k++;
            if (i >= 6)  /* N7/C8/N9 */
                kr++;
            idx[i] = n;
        } else
            idx[i] = 0;  /* make explicit assignment */
    }

    if (k >= 3) {  /* at least 3 ring atoms */
        C1_prime = find_1st_atom(" C1'", AtomName, ib, ie, "");
        rmsd = check_nt_type_by_rmsd(idx, C1_prime, xyz);
    }

    if (rmsd != DUMMY && rmsd <= Gvars.NT_CUTOFF) {
        if (kr)
            return 1;
        set_U_C5M(AtomName, xyz, ib, ie);
        return 0;
    }

    CA = find_1st_atom(" CA ", AtomName, ib, ie, "");
    C = find_1st_atom(" C  ", AtomName, ib, ie, "");
    if (!C)  /* if C does not exist, use N */
        C = find_1st_atom(" N  ", AtomName, ib, ie, "");
    if (CA && C && Miscs[CA][0] == 'A' && Miscs[C][0] == 'A' &&  /* ATOM record */
        within_limits(xyz[CA], xyz[C], 0, dcrt))
        return -1;

    return -2;
}

void normalize_atom_symbol(char *asym)
{
    upperstr(asym);

    if (strlen(asym) == 1) {  /* "H" --> " H" */
        asym[1] = asym[0];
        asym[0] = ' ';
        asym[2] = '\0';
    }

    if (is_equal_string(asym, " D"))
        strcpy(asym, " H");
}

/* get the correspondence between 4-letter atom name to atomic symbol */
void get_atomlist(char **atomlist, long *num_sa)
{
    char BDIR[BUF512], str[BUF512], aname4[BUF512], asym[BUF512];
    long n = 0;
    FILE *fp;

    get_BDIR(BDIR, ATOM_FILE);
    strcat(BDIR, ATOM_FILE);
    if (Gvars.VERBOSE)
        fprintf(stderr, " ...... reading file: %s ...... \n", ATOM_FILE);

    fp = open_file(BDIR, "r");

    while ((fgets(str, sizeof str, fp) != NULL)) {
        if (str[0] == '#')
            continue;
        if (sscanf(str, "%s %s", aname4, asym) != 2)
            continue;
        if (aname4[0] == '#' || asym[0] == '#')
            continue;

        if (strlen(aname4) != 4) {
            if (Gvars.VERBOSE)
                fprintf(stderr, "atom name must be 4-char long with only [.A-Z]: <%s>\n", aname4);
            continue;
        }

        if (strlen(asym) != 1 && strlen(asym) != 2) {
            if (Gvars.VERBOSE)
                fprintf(stderr, "atomic symbol must be 1/2-char with only [A-Z]: <%s>\n", asym);
            continue;
        }

        normalize_atom_symbol(asym);
        if (!num_strmatch(asym, Gvars.ATOM_NAMES, 0, Gvars.NUM_ELE)) {
            if (Gvars.VERBOSE)
                fprintf(stderr, "skip invalid atom symbol <%s : %s>\n", asym, aname4);
            continue;
        }

        upperstr(aname4);

        if (++n > BUFBIG)
            fatal("too many atom types\n");

        sprintf(atomlist[n], "%4.4s%2.2s", aname4, asym);
    }

    close_file(fp);

    *num_sa = n;
}

long has_atom_name(long ib, long ie, char **AtomName, char *aname)
{
    long i;

    for (i = ib; i <= ie; i++)
        if (strcmp(AtomName[i], aname) == 0)
            return 1L;

    return 0L;
}

/* get one letter base name according to 3-letter to 1-letter conversion */
static char base_ident(long ib, long ie, char **AtomName, double **xyz, long isR,
                       char *rname, char *idmsg, long num_sb, char **baselist)
{
    char c = 'X';
    long i;

    for (i = 1; i <= num_sb; i++)
        if (!strncmp(rname, baselist[i], 3))
            break;

    if (i > num_sb) {
        if (isR) {
            if (has_atom_name(ib, ie, AtomName, " O6 "))  /* instead of N2: 2013oct13/14 */
                c = 'g';
            else if (!has_atom_name(ib, ie, AtomName, " N6 ") &&
                     has_atom_name(ib, ie, AtomName, " N2 "))
                c = 'g';
            else
                c = 'a';
        } else {
            if (has_atom_name(ib, ie, AtomName, " N4 "))
                c = 'c';
            else if (has_atom_name(ib, ie, AtomName, " C5M"))
                c = 't';
            else
                c = 'u';
            {
                long c1p, n1, c5;
                c1p = find_1st_atom(" C1'", AtomName, ib, ie, "");
                n1 = find_1st_atom(" N1 ", AtomName, ib, ie, "");
                c5 = find_1st_atom(" C5 ", AtomName, ib, ie, "");
                if (c1p && n1 && c5 && p1p2_dist(xyz[c1p], xyz[n1]) > 2.0 &&
                    p1p2_dist(xyz[c1p], xyz[c5]) <= 2.0)
                    c = 'p';  /* as in 3TD of 5afi */
            }
        }
        fprintf(stderr, "Match '%s' to '%c' for %s\n", rname, c, idmsg);
        fprintf(stderr, "    check it & consider to add line '%s     %c' to file <%s>\n",
                rname, c, BASE_FILE);
    } else {
        c = baselist[i][3];
        if (!isupper((int) c) || (c == 'P'))
            fprintf(stderr, "[i] uncommon %s assigned to: %c\n", idmsg, c);
    }

    return c;
}

/* get the correspondence between 3-letter and 1-letter base names */
void get_baselist(char **baselist, long *num_sb)
{
    char BDIR[BUF512], str[BUF512], base3[BUF512], base1[BUF512];
    long n = 0;
    FILE *fp;

    get_BDIR(BDIR, BASE_FILE);
    strcat(BDIR, BASE_FILE);
    if (Gvars.VERBOSE)
        fprintf(stderr, " ...... reading file: %s ...... \n", BASE_FILE);

    fp = open_file(BDIR, "r");

    while ((fgets(str, sizeof str, fp) != NULL)) {
        if (str[0] == '#')
            continue;
        if (sscanf(str, "%s %s", base3, base1) != 2)
            continue;
        if (base3[0] == '#' || base1[0] == '#')
            continue;

        if (strlen(base3) > 3 || strlen(base1) != 1) {
            if (Gvars.VERBOSE)
                fprintf(stderr, "ignoring unacceptable format: %s\n", str);
            continue;
        }

        if (++n > BUFBIG)
            fatal("too many base types\n");

        upperstr(base3);
        sprintf(baselist[n], "%3.3s%c", base3, base1[0]);
    }

    close_file(fp);

    *num_sb = n;
}

/* get base sequence, similar to GET_BPSEQ */
void get_seq(long num_residue, long **seidx, char **AtomName, char **ResName,
             char *ChainID, long *ResSeq, char **Miscs, double **xyz, char *bseq, long *RY)
{
    char idmsg[BUF512];
    long i, ib, ie;

    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];
        RY[i] = residue_ident(AtomName, xyz, Miscs, ib, ie);
        if (RY[i] >= 0) {
            sprintf(idmsg, "residue %3s %4ld%c on chain %c [#%ld]",
                    ResName[ib], ResSeq[ib], Miscs[ib][2], ChainID[ib], i);
            bseq[i] = base_ident(ib, ie, AtomName, xyz, RY[i], ResName[ib], idmsg,
                                 Gvars.NUM_SBASE, Gvars.BASELIST);
        }
    }
}

/* get base (pair) sequence, change (+) modified residue to lower case */
void get_bpseq(long ds, long num_bp, long **pair_num, long **seidx, char **AtomName,
               char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
               char **bp_seq, long *RY)
{
    char idmsg[BUF512];
    long i, ib, ie, j, rnum;

    for (i = 1; i <= ds; i++) {
        for (j = 1; j <= num_bp; j++) {
            rnum = pair_num[i][j];
            ib = seidx[rnum][1];
            ie = seidx[rnum][2];
            RY[rnum] = residue_ident(AtomName, xyz, Miscs, ib, ie);
            sprintf(idmsg, "residue %3s %4ld%c on chain %c [#%ld]",
                    ResName[ib], ResSeq[ib], Miscs[ib][2], ChainID[ib], rnum);
            if (RY[rnum] >= 0)
                bp_seq[i][j] = base_ident(ib, ie, AtomName, xyz, RY[rnum], ResName[ib],
                                          idmsg, Gvars.NUM_SBASE, Gvars.BASELIST);
            else
                fatal("Non-base: %s\n", idmsg);
        }
    }
}

/* find the 1st matching index, or return NO_MATCH */
long strmatch_idx(char *str, char **strmat, long nb, long ne)
{
    long i;

    for (i = nb; i <= ne; i++)
        if (!strcmp(str, strmat[i]))
            return i;

    return NO_MATCH;
}

/*  return number of matchs of str in strmat */
long num_strmatch(char *str, char **strmat, long nb, long ne)
{
    long i, num = 0;

    for (i = nb; i <= ne; i++)
        if (!strcmp(str, strmat[i]))
            num++;

    return num;
}

void get_idmsg(char *rname, char cid, long snum, char icode, char *idmsg)
{
    sprintf(idmsg, ": residue name '%s', chain %c, number [%4ld%c]", rname, cid, snum, icode);
}

/* return index of the first match, or 0 for no-match */
long find_1st_atom(char *str, char **strmat, long nb, long ne, char *idmsg)
{
    char aname[BUF32];
    long i, num;

    num = num_strmatch(str, strmat, nb, ne);

    if (!num) {
        if (strcmp(idmsg, "")) {
            strcpy(aname, str);
            if (Gvars.PDBV3) {
                if (is_equal_string(" O1P", str))
                    strcpy(aname, " OP1");
                else if (is_equal_string(" O2P", str))
                    strcpy(aname, " OP2");
                else if (is_equal_string(" C5M", str))
                    strcpy(aname, " C7 ");
            }
            fprintf(stderr, "[i] missing '%s' atom %s\n", aname, idmsg);
        }
        return 0;
    }

    if (num > 1 && strcmp(idmsg, "")) {
        fprintf(stderr, "more than one %s atoms %s\n", str, idmsg);
        fprintf(stderr, "   *****the first atom is used*****\n");
    }

    for (i = nb; i <= ne; i++)
        if (!strcmp(str, strmat[i]))
            break;

    return i;
}

/* get torsion angle a-b-c-d in degrees */
double torsion(double **d)
{
    double ang_deg, **vec3;
    long i;

    vec3 = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= 3; i++) {
        ddxyz(d[i], d[i + 1], vec3[i]);
        if (veclen(vec3[i]) > BOND_UPPER_LIMIT) {
            ang_deg = EMPTY_NUMBER;  /* not directly linked */
            goto RTN_HERE;
        }
    }
    negate_xyz(vec3[1]);  /* b-->a */
    ang_deg = vec_ang(vec3[1], vec3[3], vec3[2]);

  RTN_HERE:
    free_dmatrix(vec3, 1, 3, 1, 3);

    return ang_deg;
}

/* get torsion angle a-b-c-d in degrees: no consideration of breaks */
double torsion2(double **d)
{
    double ang_deg, **vec3;
    long i;

    vec3 = dmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++)
        ddxyz(d[i], d[i + 1], vec3[i]);
    negate_xyz(vec3[1]);  /* b-->a */
    ang_deg = vec_ang(vec3[1], vec3[3], vec3[2]);
    free_dmatrix(vec3, 1, 3, 1, 3);

    return ang_deg;
}

/* search the directory containing standard bases & parameter files:
   (1) current directory
   (2) directory defined by the environmental variable "X3DNA"
   (3) standard directory $HOME/x3dna-v2.3/config/ */
void get_BDIR(char *BDIR, char *filename)
{
    FILE *fp;

    fp = fopen(filename, "r");  /* check current directory */
    if (fp != NULL)
        strcpy(BDIR, "./");  /* current directory */
    else
        sprintf(BDIR, "%sconfig/", Gvars.X3DNA_HOMEDIR);
    close_file(fp);

    if (Gvars.VERBOSE)
        fprintf(stderr, "\n ...... %s%s...... \n", BDIR, filename);
}

/* transform xyz coordinates so that "haxis" is aligned to the z-axis */
void align2zaxis(long num, double *haxis, double **rotmat, double **xyz, double **xyzH)
{
    double z[4] = { EMPTY_NUMBER, 0.0, 0.0, 1.0 };
    double ang_deg, hinge[4], **rotmatT;

    rotmatT = dmatrix(1, 3, 1, 3);

    cross(haxis, z, hinge);
    ang_deg = magang(haxis, z);
    arb_rotation(hinge, ang_deg, rotmat);
    transpose_matrix(rotmat, 3, 3, rotmatT);
    multi_matrix(xyz, num, 3, rotmatT, 3, 3, xyzH);

    free_dmatrix(rotmatT, 1, 3, 1, 3);
}

/* calculate the covariance matrix between two matrices */
void cov_matrix(double **a, double **b, long nr, long nc, double **cmtx)
{
    double ave_a[4], ave_b[4];
    double **ta, **ta_x_b;
    long i, j;

    ave_dmatrix(a, nr, nc, ave_a);
    ave_dmatrix(b, nr, nc, ave_b);

    ta = dmatrix(1, nc, 1, nr);  /* transpose of a */
    ta_x_b = dmatrix(1, nc, 1, nc);  /* transpose-a multiply b */

    transpose_matrix(a, nr, nc, ta);
    multi_matrix(ta, nc, nr, b, nr, nc, ta_x_b);

    for (i = 1; i <= nc; i++)
        for (j = 1; j <= nc; j++)
            cmtx[i][j] = (ta_x_b[i][j] - ave_a[i] * ave_b[j] * nr) / (nr - 1);

    free_dmatrix(ta, 1, nc, 1, nr);
    free_dmatrix(ta_x_b, 1, nc, 1, nc);
}

/* least-squares fitting between two structures */
double ls_fitting(double **sxyz, double **exyz, long n, double **fitted_xyz, double **R,
                  double *orgi)
{
    double temp, rms_value;
    double ave_exyz[4], ave_sxyz[4], D[5];
    double **N, **U, **V;
    long i, j;

    if (n < 3)
        fatal("too few atoms for least-squares fitting\n");

    /* get the covariance matrix U */
    U = dmatrix(1, 3, 1, 3);
    cov_matrix(sxyz, exyz, n, 3, U);

    /* get 4-by-4 symmetric matrix N */
    N = dmatrix(1, 4, 1, 4);
    N[1][1] = U[1][1] + U[2][2] + U[3][3];
    N[2][2] = U[1][1] - U[2][2] - U[3][3];
    N[3][3] = -U[1][1] + U[2][2] - U[3][3];
    N[4][4] = -U[1][1] - U[2][2] + U[3][3];
    N[1][2] = U[2][3] - U[3][2];
    N[2][1] = N[1][2];
    N[1][3] = U[3][1] - U[1][3];
    N[3][1] = N[1][3];
    N[1][4] = U[1][2] - U[2][1];
    N[4][1] = N[1][4];
    N[2][3] = U[1][2] + U[2][1];
    N[3][2] = N[2][3];
    N[2][4] = U[3][1] + U[1][3];
    N[4][2] = N[2][4];
    N[3][4] = U[2][3] + U[3][2];
    N[4][3] = N[3][4];

    /* get N's eigenvalues and eigenvectors */
    V = dmatrix(1, 4, 1, 4);
    jacobi(N, 4, D, V);

    /* get the rotation matrix */
    for (i = 1; i <= 4; i++)
        for (j = 1; j <= 4; j++)
            N[i][j] = V[i][4] * V[j][4];
    R[1][1] = N[1][1] + N[2][2] - N[3][3] - N[4][4];
    R[1][2] = 2 * (N[2][3] - N[1][4]);
    R[1][3] = 2 * (N[2][4] + N[1][3]);
    R[2][1] = 2 * (N[3][2] + N[1][4]);
    R[2][2] = N[1][1] - N[2][2] + N[3][3] - N[4][4];
    R[2][3] = 2 * (N[3][4] - N[1][2]);
    R[3][1] = 2 * (N[4][2] - N[1][3]);
    R[3][2] = 2 * (N[4][3] + N[1][2]);
    R[3][3] = N[1][1] - N[2][2] - N[3][3] + N[4][4];

    ave_dmatrix(sxyz, n, 3, ave_sxyz);
    ave_dmatrix(exyz, n, 3, ave_exyz);

    /* fitted sxyz origin */
    for (i = 1; i <= 3; i++)
        orgi[i] = ave_exyz[i] - dot(ave_sxyz, R[i]);

    /* fitted sxyz coordinates */
    for (i = 1; i <= n; i++)
        for (j = 1; j <= 3; j++)
            fitted_xyz[i][j] = dot(sxyz[i], R[j]) + orgi[j];

    /* rms deviation */
    temp = 0.0;
    for (i = 1; i <= n; i++) {
        ddxyz(fitted_xyz[i], exyz[i], D);
        temp += dot(D, D);
    }
    rms_value = sqrt(temp / n);

    free_dmatrix(U, 1, 3, 1, 3);
    free_dmatrix(N, 1, 4, 1, 4);
    free_dmatrix(V, 1, 4, 1, 4);

    return rms_value;
}

/* fit a plane to a set of points by least squares */
void ls_plane(double **bxyz, long n, double *pnormal, double *ppos, double *odist, double *adist)
{
    double D[4];
    double **cov_mtx, **identityV, **V;
    long i, j, nml = 0;

    if (n < 3)
        fatal("too few atoms for least-squares fitting\n");

    cov_mtx = dmatrix(1, 3, 1, 3);
    V = dmatrix(1, 3, 1, 3);
    identityV = dmatrix(1, 3, 1, 3);

    cov_matrix(bxyz, bxyz, n, 3, cov_mtx);
    jacobi(cov_mtx, 3, D, V);

    identity_matrix(identityV, 3);
    for (i = 1; i <= 3 && !nml; i++)
        for (j = 1; j <= 3 && !nml; j++)
            if (fabs(V[i][j] - identityV[i][j]) > XEPS)
                nml = 1;
    if (nml)
        for (i = 1; i <= 3; i++)
            pnormal[i] = V[i][1];
    else {  /* V is an identity matrix */
        pnormal[1] = 0.0;
        pnormal[2] = 0.0;
        pnormal[3] = 1.0;
    }

    ave_dmatrix(bxyz, n, 3, ppos);

    /* make the z-component of pnormal to be positive */
    if (pnormal[3] < 0)
        negate_xyz(pnormal);

    /* distance from the origin to the plane */
    *odist = dot(ppos, pnormal);

    /* distance from each point to the plane */
    for (i = 1; i <= n; i++)
        adist[i] = dot(bxyz[i], pnormal) - *odist;

    free_dmatrix(cov_mtx, 1, 3, 1, 3);
    free_dmatrix(V, 1, 3, 1, 3);
    free_dmatrix(identityV, 1, 3, 1, 3);
}

/* get the arbitrary rotation matrix */
void arb_rotation(double *va, double ang_deg, double **rot_mtx)
{
    double c, dc, s, vlen;
    long i;

    vlen = veclen(va);
    if (vlen < XEPS)  /* [0 0 0] */
        identity_matrix(rot_mtx, 3);
    else {
        for (i = 1; i <= 3; i++)
            va[i] /= vlen;  /* unit vector */
        ang_deg = deg2rad(ang_deg);
        c = cos(ang_deg);
        s = sin(ang_deg);
        dc = 1 - c;
        rot_mtx[1][1] = c + dc * va[1] * va[1];
        rot_mtx[1][2] = va[1] * va[2] * dc - va[3] * s;
        rot_mtx[1][3] = va[1] * va[3] * dc + va[2] * s;
        rot_mtx[2][1] = va[1] * va[2] * dc + va[3] * s;
        rot_mtx[2][2] = c + dc * va[2] * va[2];
        rot_mtx[2][3] = va[2] * va[3] * dc - va[1] * s;
        rot_mtx[3][1] = va[1] * va[3] * dc - va[2] * s;
        rot_mtx[3][2] = va[2] * va[3] * dc + va[1] * s;
        rot_mtx[3][3] = c + dc * va[3] * va[3];
    }
}

/* angle in degrees between va and vb with vref for sign control
   va & vb are unchanged by making an additional copy of each
   all three vectors are 1-by-3 */
double vec_ang(double *va, double *vb, double *vref)
{
    double ang_deg, va_cp[4], vb_cp[4];

    /* make a copy of va and vb */
    cpxyz(va, va_cp);
    cpxyz(vb, vb_cp);

    /* get orthogonal components */
    vec_orth(va_cp, vref);
    vec_orth(vb_cp, vref);

    /* angle in absolute sense */
    ang_deg = magang(va_cp, vb_cp);

    if (sign_control(va_cp, vb_cp, vref) < 0)
        ang_deg = -ang_deg;

    return ang_deg;
}

/* get the vector which has certain angle with another vector */
void get_vector(double *va, double *vref, double deg_ang, double *vo)
{
    double va_cp[4], **temp;
    long i;

    cpxyz(va, va_cp);  /* make a copy of va: <ana_fncs.c> */
    if (dot(va_cp, vref) > XEPS) {
        fprintf(stderr, "Angle between va/vref: %.3f degrees\n", magang(va_cp, vref));
        vec_orth(va_cp, vref);
    }

    temp = dmatrix(1, 3, 1, 3);

    arb_rotation(vref, deg_ang, temp);
    for (i = 1; i <= 3; i++)
        vo[i] = dot(temp[i], va_cp);
    vec_norm(vo);

    free_dmatrix(temp, 1, 3, 1, 3);
}

void rotate(double **a, long i, long j, long k, long l, double *g, double *h, double s, double tau)
{
    *g = a[i][j];
    *h = a[k][l];
    a[i][j] = *g - s * (*h + *g * tau);
    a[k][l] = *h + s * (*g - *h * tau);
}

/* sort eigenvalues into ascending order and rearrange eigenvectors */
void eigsrt(double *d, double **v, long n)
{
    double p;
    long i, j, k;

    for (i = 1; i < n; i++) {
        p = d[k = i];
        for (j = i + 1; j <= n; j++)
            if (d[j] < p)
                p = d[k = j];
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 1; j <= n; j++)
                dval_swap(&v[j][i], &v[j][k]);
        }
    }
}

void jacobi(double **a, long n, double *d, double **v)
{
    long i, j, iq, ip;
    double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

    b = dvector(1, n);
    z = dvector(1, n);
    identity_matrix(v, n);
    for (ip = 1; ip <= n; ip++) {
        b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }
    for (i = 1; i <= 100; i++) {
        sm = 0.0;
        for (ip = 1; ip <= n - 1; ip++) {
            for (iq = ip + 1; iq <= n; iq++)
                sm += fabs(a[ip][iq]);
        }
        if (sm < XEPS) {
            free_dvector(z, 1, n);
            free_dvector(b, 1, n);
            eigsrt(d, v, n);
            return;
        }
        if (i < 4)
            tresh = 0.2 * sm / (n * n);
        else
            tresh = 0.0;
        for (ip = 1; ip <= n - 1; ip++) {
            for (iq = ip + 1; iq <= n; iq++) {
                g = 100.0 * fabs(a[ip][iq]);
                if (i > 4 && (fabs(d[ip]) + g) == fabs(d[ip])
                    && (fabs(d[iq]) + g) == fabs(d[iq]))
                    a[ip][iq] = 0.0;
                else if (fabs(a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if ((fabs(h) + g) == fabs(h))
                        t = a[ip][iq] / h;
                    else {
                        theta = 0.5 * h / a[ip][iq];
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0)
                            t = -t;
                    }
                    c = 1.0 / sqrt(1 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 1; j <= ip - 1; j++)
                        rotate(a, j, ip, j, iq, &g, &h, s, tau);
                    for (j = ip + 1; j <= iq - 1; j++)
                        rotate(a, ip, j, j, iq, &g, &h, s, tau);
                    for (j = iq + 1; j <= n; j++)
                        rotate(a, ip, j, iq, j, &g, &h, s, tau);
                    for (j = 1; j <= n; j++)
                        rotate(v, j, ip, j, iq, &g, &h, s, tau);
                }
            }
        }
        for (ip = 1; ip <= n; ip++) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    fatal("too many iterations\n");
}

void dludcmp(double **a, long n, long *indx, double *d)
{
    double big, dum, sum, temp;
    double *vv;
    long i, j, k;
    long imax = 0;  /* initialization */

    vv = dvector(1, n);
    *d = 1.0;
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if ((temp = fabs(a[i][j])) > big)
                big = temp;
        if (big == 0.0)
            fatal("singular matrix in routine dludcmp\n");
        vv[i] = 1.0 / big;
    }
    for (j = 1; j <= n; j++) {
        for (i = 1; i < j; i++) {
            sum = a[i][j];
            for (k = 1; k < i; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 1; k <= n; k++)
                dval_swap(&a[imax][k], &a[j][k]);
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0)
            a[j][j] = XEPS;
        if (j != n) {
            dum = 1.0 / a[j][j];
            for (i = j + 1; i <= n; i++)
                a[i][j] *= dum;
        }
    }
    free_dvector(vv, 1, n);
}

void dlubksb(double **a, long n, long *indx, double *b)
{
    double sum;
    long i, ii = 0, ip, j;

    for (i = 1; i <= n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
            for (j = ii; j <= i - 1; j++)
                sum -= a[i][j] * b[j];
        else if (sum)
            ii = i;
        b[i] = sum;
    }
    for (i = n; i >= 1; i--) {
        sum = b[i];
        for (j = i + 1; j <= n; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

void dinverse(double **a, long n, double **y)
{
    double d, *col;
    long i, j, *indx;

    col = dvector(1, n);
    indx = lvector(1, n);

    dludcmp(a, n, indx, &d);

    for (j = 1; j <= n; j++) {
        for (i = 1; i <= n; i++)
            col[i] = 0.0;
        col[j] = 1.0;
        dlubksb(a, n, indx, col);
        for (i = 1; i <= n; i++)
            y[i][j] = col[i];
    }

    free_lvector(indx, 1, n);
    free_dvector(col, 1, n);
}

void rotx(double ang_deg, double **rotmat)
{
    double ang, c, s;

    ang = deg2rad(ang_deg);
    c = cos(ang);
    s = sin(ang);

    rotmat[1][1] = 1.0;
    rotmat[1][2] = 0.0;
    rotmat[1][3] = 0.0;
    rotmat[2][1] = 0.0;
    rotmat[2][2] = c;
    rotmat[2][3] = -s;
    rotmat[3][1] = 0.0;
    rotmat[3][2] = s;
    rotmat[3][3] = c;
}

void roty(double ang_deg, double **rotmat)
{
    double ang, c, s;

    ang = deg2rad(ang_deg);
    c = cos(ang);
    s = sin(ang);

    rotmat[1][1] = c;
    rotmat[1][2] = 0.0;
    rotmat[1][3] = s;
    rotmat[2][1] = 0.0;
    rotmat[2][2] = 1.0;
    rotmat[2][3] = 0.0;
    rotmat[3][1] = -s;
    rotmat[3][2] = 0.0;
    rotmat[3][3] = c;
}

void rotz(double ang_deg, double **rotmat)
{
    double ang, c, s;

    ang = deg2rad(ang_deg);
    c = cos(ang);
    s = sin(ang);

    rotmat[1][1] = c;
    rotmat[1][2] = -s;
    rotmat[1][3] = 0.0;
    rotmat[2][1] = s;
    rotmat[2][2] = c;
    rotmat[2][3] = 0.0;
    rotmat[3][1] = 0.0;
    rotmat[3][2] = 0.0;
    rotmat[3][3] = 1.0;
}

void get_alc_nums(char *alcname, long *num, long *nbond)
{
    char str[BUF512];
    FILE *fp;

    fp = open_file(alcname, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld %*s %ld", num, nbond) != 2)
        fatal("can not read atom & bond numbers: %s\n", str);

    if (!*num)
        fprintf(stderr, "ALCHEMY file <%s> has NO atoms\n", alcname);

    close_file(fp);
}

void read_alc(char *alcname, long *num, long *nbond, char **AtomName, double **xyz,
              long *ibase, long **linkage)
{
    char str[BUF512], temp[BUF512];
    double c;
    long i;
    FILE *fp;

    fp = open_file(alcname, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld %*s %ld", num, nbond) != 2)
        fatal("cannot read atom & bond numbers: %s\n", str);

    for (i = 1; i <= *num; i++)
        if (fgets(str, sizeof str, fp) != NULL) {
            strncpy(AtomName[i], str + 6, 2);
            AtomName[i][2] = '\0';
            upperstr(AtomName[i]);  /* to upper case */
            strncpy(temp, str + 11, 9);
            temp[9] = '\0';
            if (sscanf(temp, "%lf", &xyz[i][1]) != 1)
                fatal("error reading x-coordinate\n");
            strncpy(temp, str + 20, 9);
            temp[9] = '\0';
            if (sscanf(temp, "%lf", &xyz[i][2]) != 1)
                fatal("error reading y-coordinate\n");
            strncpy(temp, str + 29, 9);
            temp[9] = '\0';
            if (sscanf(temp, "%lf", &xyz[i][3]) != 1)
                fatal("error reading z-coordinate\n");
            strncpy(temp, str + 40, 9);
            temp[9] = '\0';
            if (sscanf(temp, "%lf", &c) == 1)
                ibase[i] = lround(10.0 * c);
            else
                ibase[i] = NON_WC_IDX;  /* default */
        } else
            fatal("error in reading atom records\n");

    for (i = 1; i <= *nbond; i++)
        if (fgets(str, sizeof str, fp) != NULL) {
            strncpy(temp, str + 6, 5);
            temp[5] = '\0';
            if (sscanf(temp, "%ld", &linkage[i][1]) != 1)
                fatal("error reading linkage atom 1\n");
            strncpy(temp, str + 12, 5);
            temp[5] = '\0';
            if (sscanf(temp, "%ld", &linkage[i][2]) != 1)
                fatal("error reading linkage atom 2\n");
        } else
            fatal("error in reading linkage information\n");

    close_file(fp);
}

void write_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
               long **linkage, char *alcfile)
{
    char aname[5];
    long i;
    FILE *fp;

    reset_xyz(num, xyz, "f9.4");

    fp = open_file(alcfile, "w");
    fprintf(fp, "%5ld ATOMS, %5ld BONDS\n", num, nbond);
    for (i = 1; i <= num; i++) {
        if (sscanf(AtomName[i], "%s", aname) == EOF)
            strcpy(aname, "C ");
        if (isdigit((int) aname[0])) {  /* like "5H1 " etc */
            aname[0] = aname[1];
            aname[1] = aname[2];
        }
        if (isdigit((int) aname[1]) || aname[1] == '\0')
            aname[1] = ' ';
        aname[2] = '\0';
        fprintf(fp, "%5ld %-2s   %9.4f%9.4f%9.4f  %9.4f\n", i, aname, xyz[i][1],
                xyz[i][2], xyz[i][3], (ibase == NULL) ? 0.0 : ibase[i] / 10.0);
    }
    for (i = 1; i <= nbond; i++)
        fprintf(fp, "%5ld %5ld %5ld\n", i, linkage[i][1], linkage[i][2]);
    close_file(fp);
}

void free_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
              long zero_1, long **linkage)
{
    free_cmatrix(AtomName, 1, num, 0, 2);
    free_dmatrix(xyz, 1, num, 1, 3);
    free_lvector(ibase, 1, num);
    free_lmatrix(linkage, 1, nbond, zero_1, 2);
}

/* add lines connecting consecutive origins to ALCHEMY file */
void cnct_org(long num_bp, long ia, long ib, char **tAtomName, double **txyz,
              long *tibase, long **tlinkage, double **org_xyz)
{
    long i, ik;

    for (i = 1; i <= num_bp; i++) {
        ik = ia + i;
        strcpy(tAtomName[ik], "O");
        cpxyz(org_xyz[i], txyz[ik]);
        tibase[ik] = NON_WC_IDX;  /* non-common base atoms */
    }

    for (i = 1; i < num_bp; i++) {
        ik = ib + i;
        tlinkage[ik][1] = ia + i;
        tlinkage[ik][2] = ia + i + 1;
    }
}

/* sort a double vector into ascending order & keep the index
   use Shell's method as in NR in C book Ed. 2, pp.331-332
   cf. lsort below */
void dsort(long n, double *a, long *idx)
{
    double v;
    long i, inc, iv, j;

    inc = 1;
    do {
        inc *= 3;
        inc++;
    } while (inc <= n);

    for (i = 1; i <= n; i++)
        idx[i] = i;

    do {
        inc /= 3;
        for (i = inc + 1; i <= n; i++) {
            v = a[i];
            iv = idx[i];
            j = i;
            while (a[j - inc] > v) {
                a[j] = a[j - inc];
                idx[j] = idx[j - inc];
                j -= inc;
                if (j <= inc)
                    break;
            }
            a[j] = v;
            idx[j] = iv;
        }
    } while (inc > 1);
}

/* sort a long vector into ascending order & keep the index
   use Shell's method as in NR in C book Ed. 2, pp.331-332 */
void lsort(long n, long *a, long *idx)
{
    long v;
    long i, inc, iv, j;

    inc = 1;
    do {
        inc *= 3;
        inc++;
    } while (inc <= n);

    for (i = 1; i <= n; i++)
        idx[i] = i;

    do {
        inc /= 3;
        for (i = inc + 1; i <= n; i++) {
            v = a[i];
            iv = idx[i];
            j = i;
            while (a[j - inc] > v) {
                a[j] = a[j - inc];
                idx[j] = idx[j - inc];
                j -= inc;
                if (j <= inc)
                    break;
            }
            a[j] = v;
            idx[j] = iv;
        }
    } while (inc > 1);
}

/* reverse a long vector: output replaces the original */
void lreverse(long ia, long n, long *lvec)
{
    long i, *ltmp;

    ltmp = lvector(1, n);
    for (i = 1; i <= n; i++)
        ltmp[i] = lvec[n + ia - i];
    for (i = 1; i <= n; i++)
        lvec[ia + i - 1] = ltmp[i];

    free_lvector(ltmp, 1, n);
}

/* write a title to XFIG 3.2 format */
void fig_title(FILE * fp)
{
    fprintf(fp,
            "#FIG 3.2  Creator: %s\n"
            "Portrait\n"
            "Flush Left\n"
            "Inches\n" "Letter\n" "100.00\n" "Single\n" "-2\n" "1200 2\n\n", Gvars.X3DNA_VER);
}

void ps_title_cmds(FILE * fp, char *imgfile, long *bbox)
{
    char BDIR[BUF512], str[BUF512];
    char *ps_image_par = "ps_image.par";
    long i;
    time_t run_time;
    FILE *fpp;

    run_time = time(NULL);

    fprintf(fp, "%%!PS-Adobe-3.0\n");
    fprintf(fp, "%%%%Title: (%s)\n", imgfile);
    fprintf(fp, "%%%%Creator: (%s)\n", Gvars.X3DNA_VER);
    strcpy(str, ctime(&run_time));
    str[strlen(str) - 1] = '\0';  /* change \n to \0 */
    fprintf(fp, "%%%%CreationDate: (%s)\n", str);
    fprintf(fp, "%%%%Orientation: Portrait\n");
    fprintf(fp, "%%%%BoundingBox: ");
    for (i = 1; i <= 4; i++)
        fprintf(fp, "%6ld", bbox[i]);
    fprintf(fp, "\n\n");

    /* read in color parameter file */
    get_BDIR(BDIR, ps_image_par);
    strcat(BDIR, ps_image_par);
    fprintf(stderr, " ...... reading file: %s ...... \n", ps_image_par);

    fpp = open_file(BDIR, "r");
    while (fgets(str, sizeof str, fpp) != NULL)
        fprintf(fp, "%s", str);
    close_file(fpp);
}

/* reset x/y coordinates to XFIG units */
void get_fig_xy(long num, double **xyz, long nO, double **oxyz, long *urxy,
                long frame_box, FILE * fp)
{
    char *format = "%6ld%6ld";
    double paper_size[2] = { 8.5, 11.0 };  /* US letter size */
    double max_xy[3];
    long i, j, llxy[3];

    /* change y-axis according to XFIG coordinate system */
    for (i = 1; i <= num; i++)
        xyz[i][2] = urxy[2] - xyz[i][2];
    if (nO)
        for (i = 1; i <= nO; i++)
            oxyz[i][2] = urxy[2] - oxyz[i][2];

    max_dmatrix(xyz, num, 2, max_xy);
    for (i = 1; i <= 2; i++)
        urxy[i] = lround(max_xy[i]);  /* new urxy */

    /* 1/4 offset the figure on a US letter (8.5in-by-11in) */
    for (i = 1; i <= 2; i++)
        llxy[i] = lround(0.25 * (paper_size[i - 1] * 1200 - urxy[i]));
    for (i = 1; i <= num; i++)
        for (j = 1; j <= 2; j++)
            xyz[i][j] += llxy[j];
    if (nO)
        for (i = 1; i <= nO; i++)
            for (j = 1; j <= 2; j++)
                oxyz[i][j] += llxy[j];

    max_dmatrix(xyz, num, 2, max_xy);
    for (i = 1; i <= 2; i++)
        urxy[i] = lround(max_xy[i]);

    fig_title(fp);

    if (frame_box) {
        /* draw a box around the figure at the deepest depth (999) */
        fprintf(fp, "# draw a boundary box here\n");
        fprintf(fp, "2 3 0 1 0 0 999 0 -1 0.0 2 1 0 0 0 5\n");
        fprintf(fp, format, llxy[1] - FIG_BOUND, llxy[2] - FIG_BOUND);
        fprintf(fp, format, urxy[1] + FIG_BOUND, llxy[2] - FIG_BOUND);
        fprintf(fp, format, urxy[1] + FIG_BOUND, urxy[2] + FIG_BOUND);
        fprintf(fp, format, llxy[1] - FIG_BOUND, urxy[2] + FIG_BOUND);
        fprintf(fp, format, llxy[1] - FIG_BOUND, llxy[2] - FIG_BOUND);
        fprintf(fp, "\n\n");
    }
}

/* get the point position of x-, y-, or z-axis */
void get_pxy(double *xy1, double *xy2, double r, double *px, double *py)
{
    double dx, dy, dd;

    dx = xy2[1] - xy1[1];
    dy = xy2[2] - xy1[2];
    dd = sqrt(dx * dx + dy * dy);

    if (dd > XEPS) {
        *px = xy2[1] + r * dx / dd;
        *py = xy2[2] + r * dy / dd;
    } else {
        *px = xy2[1];
        *py = xy2[2];
    }
}

void alc2fig(long nobj, long *idx, long *depth, long **allobj, double **blkxyz,
             double **oxyz, long *ibase, long faces[][5], long *opts, FILE * fp)
{
    char *format = "%6.0f%6.0f";
    double dot_sep, msat, Msat, px, py;
    long is_color, same_faces, updown, mgroove;
    long dlcol, dwidth, w1, w2, line_width, join_style, cap_style, o_sides, mfcol;
    long bcol_code, i, ib, ie, ioffset8, ip, j, k, Mfill, mfill;
    long **bc_idx;  /* [2][7] -- color vs black/white, with 7 types total */

    is_color = opts[2];  /* decomposed for clarity */
    same_faces = opts[7];
    updown = opts[9];
    mgroove = opts[10];

    bc_idx = lmatrix(0, 1, 0, 6);  /* [2][7] */
    get_fig_pars(&dot_sep, &dlcol, &dwidth, &w1, &w2, bc_idx, &msat, &Msat,
                 &o_sides, &line_width, &join_style, &cap_style, &mfcol);

    /* render each object from inside to outside */
    for (i = 1; i <= nobj; i++) {
        j = idx[i];  /* object number */
        if (allobj[1][j] > 0) {  /* blocks */
            bcol_code = bc_idx[is_color][ibase[allobj[1][j]]];
            k = lround(msat * 20);
            if (mgroove)  /* could still be in color */
                mfill = 20 - k;
            else
                mfill = (bcol_code) ? 20 + k : 20 - k;
            k = lround(Msat * 20);
            if (mgroove)  /* could still be in color */
                Mfill = k;
            else
                Mfill = (bcol_code) ? 40 - k : k;

            if (!same_faces) {
                if ((!updown && !allobj[2][j]) ||  /* minor groove side */
                    (updown && allobj[2][j] == 4))  /* upper face */
                    fprintf(fp, "2 3 0 %2ld %2ld %2ld %4ld 0 %2ld",
                            line_width, bcol_code, bcol_code, depth[i], mfill);
                else if ((!updown && allobj[2][j] == 1) ||  /* major groove side: back */
                         (updown && allobj[2][j] == 5))  /* lower face */
                    fprintf(fp, "2 3 0 %2ld %2ld %2ld %4ld 0 %2ld",
                            line_width, bcol_code, bcol_code, depth[i], Mfill);
                else {  /* other sides */
                    if (bcol_code)  /* color */
                        fprintf(fp, "2 3 0 %2ld %2ld %2ld %4ld 0 20",
                                line_width, bcol_code, o_sides, depth[i]);
                    else  /* black & white */
                        fprintf(fp, "2 3 0 %2ld 0 0 %4ld 0 0", line_width, depth[i]);
                }
            } else {
                if (mgroove && !allobj[2][j])  /* minor groove side */
                    fprintf(fp, "2 3 0 %2ld 0 %2ld %4ld 0 %2ld", line_width,
                            bc_idx[0][0], depth[i], mfcol);
                else
                    fprintf(fp, "2 3 0 %2ld 0 %2ld %4ld 0 %2ld", line_width, bcol_code,
                            depth[i], (bcol_code) ? mfill : Mfill);
            }

            fprintf(fp, " 0.0 %2ld %2ld 0 0 0 5\n", join_style, cap_style);
            ioffset8 = (allobj[1][j] - 1) * 8;
            for (k = 0; k < 5; k++) {
                ip = ioffset8 + faces[allobj[2][j]][k];
                fprintf(fp, format, blkxyz[ip][1], blkxyz[ip][2]);
            }
            fprintf(fp, "\n");
        } else {  /* O connection lines */
            if (!is_color)
                dlcol = 0;  /* black */
            if (!allobj[1][j]) {  /* origin lines */
                fprintf(fp, "2 1 2 %2ld %2ld 0 %4ld 0 -1 %5.1f", dwidth, dlcol, depth[i], dot_sep);
                fprintf(fp, " %2ld %2ld 0 0 0 2\n", join_style, cap_style);
            } else if (allobj[1][j] == -1) {  /* helix axis */
                fprintf(fp, "2 1 0 %2ld %2ld 0 %4ld 0 -1 0.0", w2, dlcol, depth[i]);
                fprintf(fp, " %2ld %2ld 0 0 0 2\n", join_style, cap_style);
            } else {  /* reference axis */
                fprintf(fp, "2 1 0 %2ld %2ld 0 %4ld 0 -1 0.0", w1, dlcol, depth[i]);
                fprintf(fp, " %2ld %2ld 0 1 0 2\n", join_style, cap_style);
                fprintf(fp, "  2 0 2 36 66\n");
            }
            ib = allobj[2][j] / 10000;
            fprintf(fp, format, oxyz[ib][1], oxyz[ib][2]);
            ie = allobj[2][j] % 10000;
            fprintf(fp, format, oxyz[ie][1], oxyz[ie][2]);
            fprintf(fp, "\n");
            if (allobj[1][j] <= -2) {
                get_pxy(oxyz[ib], oxyz[ie], 10.0, &px, &py);
                fprintf(fp, "4 1 %2ld %4ld 0 18 18 0.0 4 165 165", dlcol, depth[i]);
                /* Helvetica-Bold: #18 */
                fprintf(fp, format, px, py);
                if (allobj[1][j] == -2)
                    fprintf(fp, " x\\001\n");
                else if (allobj[1][j] == -3)
                    fprintf(fp, " y\\001\n");
                else
                    fprintf(fp, " z\\001\n");
            }
        }
    }

    free_lmatrix(bc_idx, 0, 1, 0, 6);
}

/* reset x/y coordinates to PS units */
void get_ps_xy(char *imgfile, long *urxy, long frame_box, FILE * fp)
{
    char *format = "%6ld%6ld";
    double paper_size[2] = { 8.5, 11.0 };  /* US letter size */
    long i;
    long bbox[5], llxy[3];

    /* centralize the figure on a US letter (8.5in-by-11in) */
    for (i = 1; i <= 2; i++)
        llxy[i] = lround(0.5 * (paper_size[i - 1] * 72 - urxy[i]));

    /* boundary box */
    for (i = 1; i <= 2; i++) {
        bbox[i] = llxy[i] - PS_BOUND;
        bbox[i + 2] = urxy[i] + llxy[i] + PS_BOUND;
    }

    ps_title_cmds(fp, imgfile, bbox);

    fprintf(fp, "%6ld%6ld translate\n\n", llxy[1], llxy[2]);

    if (frame_box) {
        /* draw a box around the figure */
        fprintf(fp, "NP ");
        fprintf(fp, format, -PS_BOUND, -PS_BOUND);
        fprintf(fp, format, urxy[1] + PS_BOUND, -PS_BOUND);
        fprintf(fp, format, urxy[1] + PS_BOUND, urxy[2] + PS_BOUND);
        fprintf(fp, format, -PS_BOUND, urxy[2] + PS_BOUND);
        fprintf(fp, " DB stroke\n\n");
    }
}

void alc2ps(long nobj, long *idx, long **allobj, double **blkxyz, double **oxyz,
            long *ibase, long faces[][5], long *opts, FILE * fp)
{
    char bname;
    char *format = "%7.1f%7.1f";
    char bc_idx[2][8] = { "XXXXXXX",  /* black & white style: [8] to allow for '\0' */
        CX_LIST
    };  /* color-coded for ACGITUX */
    double px, py;
    long is_color, same_faces, updown, mgroove;
    long i, ib, ie, ioffset8, ip, j, k;

    is_color = opts[2];  /* decomposed for clarity */
    same_faces = opts[7];
    updown = opts[9];
    mgroove = opts[10];

    /* render each object from inside to outside */
    for (i = 1; i <= nobj; i++) {
        j = idx[i];  /* object number */
        fprintf(fp, "NP ");
        if (allobj[1][j] > 0) {  /* blocks */
            bname = bc_idx[is_color][ibase[allobj[1][j]]];
            fprintf(fp, "%cl ", (!same_faces) ? bname : 'X');
            ioffset8 = (allobj[1][j] - 1) * 8;
            for (k = 0; k < 4; k++) {
                ip = ioffset8 + faces[allobj[2][j]][k];
                fprintf(fp, format, blkxyz[ip][1], blkxyz[ip][2]);
            }
            fprintf(fp, " DB\n");
            if (!same_faces) {
                if ((!updown && !allobj[2][j]) ||  /* minor groove side */
                    (updown && allobj[2][j] == 4))  /* upper face */
                    fprintf(fp, "  gsave %cm grestore stroke\n", bname);
                else if ((!updown && allobj[2][j] == 1) ||  /* major groove side: back */
                         (updown && allobj[2][j] == 5))  /* lower face */
                    fprintf(fp, "  gsave %cM grestore stroke\n", bname);
                else {  /* other sides */
                    if (bname == 'X')  /* black & white */
                        fprintf(fp, "  gsave 1.0 FB grestore stroke\n");
                    else
                        fprintf(fp, "  gsave OTHER_SIDES grestore stroke\n");
                }
            } else {
                if (mgroove && !allobj[2][j])  /* minor groove side */
                    fprintf(fp, "  gsave Sm grestore stroke\n");
                else {
                    if (bname == 'X')  /* black & white */
                        fprintf(fp, "  gsave XM grestore stroke\n");
                    else
                        fprintf(fp, "  gsave %cm grestore stroke\n", bname);
                }
            }
        } else {  /* O connection lines */
            (is_color) ? fprintf(fp, "Dl ") : fprintf(fp, "Xl ");
            ib = allobj[2][j] / 10000;
            fprintf(fp, format, oxyz[ib][1], oxyz[ib][2]);
            ie = allobj[2][j] % 10000;
            fprintf(fp, format, oxyz[ie][1], oxyz[ie][2]);
            if (!allobj[1][j])  /* origin lines */
                fprintf(fp, "  gsave Dw Ds LN grestore\n");
            else if (allobj[1][j] == -1)  /* helix axis */
                fprintf(fp, "  gsave W2 LN grestore\n");
            else {  /* reference axis */
                fprintf(fp, "  gsave W1 LN grestore\n");

                get_pxy(oxyz[ib], oxyz[ie], 8.0, &px, &py);
                fprintf(fp, format, px, py);
                fprintf(fp, "  moveto");
                if (allobj[1][j] == -2)
                    fprintf(fp, " (x) SCENTER\n");
                else if (allobj[1][j] == -3)
                    fprintf(fp, " (y) SCENTER\n");
                else
                    fprintf(fp, " (z) SCENTER\n");
            }
        }
    }
    fprintf(fp, "\nshowpage\n");
}

/* get base ring atom index in one residue */
void bring_atoms(long ib, long ie, long ra_num, char **AtomName, long *nmatch, long *batom)
{
    static char *RingAtom[] = { RA_LIST };
    long i, j;

    *nmatch = 0;

    for (i = 0; i < ra_num; i++) {
        j = find_1st_atom(RingAtom[i], AtomName, ib, ie, "in base ring atoms");
        if (j)
            batom[++*nmatch] = j;
    }
}

/* get base ring atom index for all residues: num_ring */
void all_bring_atoms(long num_residue, long *RY, long **seidx, char **AtomName,
                     long *num_ring, long **ring_atom)
{
    long i, j, nmatch;

    *num_ring = 0;
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] < 0) {  /* non-base residue */
            ring_atom[i][10] = -1;
            continue;
        }
        j = (RY[i] == 1) ? 9 : 6;
        bring_atoms(seidx[i][1], seidx[i][2], j, AtomName, &nmatch, ring_atom[i]);
        if (nmatch == j) {
            ring_atom[i][10] = j;
            ++*num_ring;
        }
    }
}

/* get base index for coloring purpose */
void base_idx(long num, char *bseq, long *ibase, long single)
{
    static char *cmn_base = CB_LIST, *pchar;
    long i;

    if (single) {  /* for a single case */
        if ((pchar = strchr(cmn_base, toupper((int) *bseq))) != NULL)
            *ibase = pchar - cmn_base;
        else
            *ibase = NON_WC_IDX;
    } else {
        for (i = 1; i <= num; i++)
            if ((pchar = strchr(cmn_base, toupper((int) bseq[i]))) != NULL)
                ibase[i] = pchar - cmn_base;
            else
                ibase[i] = NON_WC_IDX;
    }
}

/* get basepair index for coloring purpose */
long basepair_idx(char *bpi)
{
    static char *WC[9] = { WC_LIST };
    long i, bidx;

    i = find_1st_atom(bpi, WC, 1, 8, "");
    if (!i)
        bidx = NON_WC_IDX;  /* non A-T/G-C pair */
    else if (i <= 2)  /* AT or AU */
        bidx = 0;  /* 0 for A */
    else if (i <= 4)  /* TA or UA */
        bidx = 4;  /* 4 for T */
    else if (i <= 6)  /* GC or IC */
        bidx = 2;  /* 2 for G */
    else  /* CG or CI */
        bidx = 1;  /* 1 for C */

    return bidx;
}

/* given plane normal and its center, project all coordinates onto it */
void plane_xyz(long num, double **xyz, double *ppos, double *nml, double **nxyz)
{
    long i, j;
    double temp, d[4];

    for (i = 1; i <= num; i++) {
        ddxyz(ppos, xyz[i], d);
        temp = dot(d, nml);
        for (j = 1; j <= 3; j++)
            nxyz[i][j] = ppos[j] + d[j] - temp * nml[j];
    }
}

/* project base atoms onto its least-squares plane defined by ring atoms */
void prj2plane(long num, long ra_num, char **AtomName, double **xyz, double z0, double **nxyz)
{
    double ang, temp;
    double zaxis[4] = { EMPTY_NUMBER, 0.0, 0.0, 1.0 };
    double adist[10], hinge[4], ppos[4], z[4];
    double **bxyz, **rmtx;

    long i, j, nmatch;
    long batom[10];

    /* find base ring atoms */
    bring_atoms(1, num, ra_num, AtomName, &nmatch, batom);

    /* base least-squares plane */
    bxyz = dmatrix(1, nmatch, 1, 3);

    for (i = 1; i <= nmatch; i++)
        cpxyz(xyz[batom[i]], bxyz[i]);
    ls_plane(bxyz, nmatch, z, ppos, &temp, adist);

    /* get the new set of coordinates */
    plane_xyz(num, xyz, ppos, z, nxyz);

    /* reorient the structure to make z-coordinate zero */
    rmtx = dmatrix(1, 3, 1, 3);
    if (z0) {
        cross(z, zaxis, hinge);
        ang = magang(z, zaxis);
        arb_rotation(hinge, ang, rmtx);
        for (i = 1; i <= num; i++)
            for (j = 1; j <= 3; j++)
                nxyz[i][j] = dot(xyz[i], rmtx[j]);
        for (i = 1; i <= num; i++)
            nxyz[i][3] -= nxyz[1][3];
    }
    free_dmatrix(bxyz, 1, nmatch, 1, 3);
    free_dmatrix(rmtx, 1, 3, 1, 3);
}

/* reset x & y coordinates to fit the scale */
void adjust_xy(long num, double **xyz, long nO, double **oxyz, double scale_factor,
               long default_size, long *urxy)
{
    long i, j;
    double temp;
    double max_xy[3], min_xy[3];

    max_dmatrix(xyz, num, 2, max_xy);
    min_dmatrix(xyz, num, 2, min_xy);

    /* get maximum dx or dy */
    temp = dval_max(max_xy[1] - min_xy[1], max_xy[2] - min_xy[2]);

    scale_factor = fabs(scale_factor);
    if (scale_factor < XEPS)
        scale_factor = default_size / temp;

    fprintf(stderr, "\n ...... scale factor: %.2f ...... \n", scale_factor);

    move_position(xyz, num, 2, min_xy);
    for (i = 1; i <= num; i++)
        for (j = 1; j <= 2; j++)
            xyz[i][j] *= scale_factor;
    if (nO) {
        move_position(oxyz, nO, 2, min_xy);
        for (i = 1; i <= nO; i++)
            for (j = 1; j <= 2; j++)
                oxyz[i][j] *= scale_factor;
    }
    max_dmatrix(xyz, num, 2, max_xy);
    for (i = 1; i <= 2; i++)
        urxy[i] = lround(max_xy[i]);
}

void get_depth(long nobj, long *zval, long *depth)
{
    double temp;
    long depth_low = 991, depth_up = 11;  /* depth level */
    long i, j;

    /* reset zval to [depth_low -- depth_up] for depth level */
    j = depth_low - depth_up;
    temp = zval[nobj] - zval[1];
    for (i = 1; i <= nobj; i++)
        depth[i] = lround(depth_low - j * (zval[i] - zval[1]) / temp);
}

void raster3d_header(long num, double **xyz, double scale_factor, long no_header,
                     long frame_box, FILE * fp)
{
    char BDIR[BUF512], str[BUF512], *header_file = "my_header.r3d";
    double temp, ave_xyz[4], min_xyz[4], max_xyz[4];
    double rad = 0.06, rgbv[4] = { 0.0, 0.25, 0.25, 0.25 };
    long i, itype = 3;
    FILE *fpp;

    min_dmatrix(xyz, num, 3, min_xyz);
    max_dmatrix(xyz, num, 3, max_xyz);
    avexyz(max_xyz, min_xyz, ave_xyz);
    temp = dval_max(max_xyz[1] - min_xyz[1], max_xyz[2] - min_xyz[2]);

    if (!no_header) {  /* write header section */
        if (scale_factor < XEPS)
            scale_factor = 6.0 + temp;
        fprintf(stderr, "\n ...... scale factor: %.2f ...... \n", scale_factor);

        get_BDIR(BDIR, header_file);
        strcat(BDIR, header_file);
        fprintf(stderr, " ...... reading file: %s ...... \n", header_file);

        fpp = open_file(BDIR, "r");
        for (i = 1; i <= 20; i++) {
            if (fgets(str, sizeof str, fpp) == NULL)
                fatal("error reading header.r3d\n");
            (i != 16) ? fputs(str, fp) :
                fprintf(fp, "%9.3f%9.3f%9.3f%9.3f\n", -ave_xyz[1],
                        -ave_xyz[2], -ave_xyz[3], scale_factor);
        }
        close_file(fpp);
    }

    if (frame_box) {
        fprintf(fp, "### Section of the frame box: 4 lines\n");
        ave_xyz[3] = min_xyz[3] = max_xyz[3];
        /* box strongly influenced by z-coordinate */
        ave_xyz[1] = max_xyz[1];  /* lower-right corner */
        ave_xyz[2] = min_xyz[2];
        r3d_rod(itype, min_xyz, ave_xyz, rad, rgbv, fp);
        r3d_rod(itype, ave_xyz, max_xyz, rad, rgbv, fp);
        ave_xyz[1] = min_xyz[1];  /* upper-left corner */
        ave_xyz[2] = max_xyz[2];
        r3d_rod(itype, min_xyz, ave_xyz, rad, rgbv, fp);
        r3d_rod(itype, ave_xyz, max_xyz, rad, rgbv, fp);
    }
}

/* read in parameters for Raster3D input */
void get_r3dpars(double **base_col, double *hb_col, double *width3, double **atom_col,
                 char *label_style)
{
    char BDIR[BUF512], str[BUF512], *raster3d_par = "raster3d.par";
    char *format = "%lf %lf %lf", *format4 = "%lf %lf %lf %lf";
    long i;
    FILE *fp;

    get_BDIR(BDIR, raster3d_par);
    strcat(BDIR, raster3d_par);
    fp = open_file(BDIR, "r");
    if (Gvars.VERBOSE)
        fprintf(stderr, " ...... reading file: %s ...... \n", raster3d_par);

    if (fgets(str, sizeof str, fp) == NULL)  /* skip one line */
        fatal("error in reading comment line\n");
    for (i = 0; i <= NBASECOL; i++)
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &base_col[i][1], &base_col[i][2], &base_col[i][3]) != 3)
            fatal("error reading base residue RGB color\n");
    if (fgets(str, sizeof str, fp) == NULL)  /* skip one line */
        fatal("error in reading comment line\n");
    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, format4, &hb_col[1], &hb_col[2], &hb_col[3], &hb_col[4]) != 4)
        fatal("error reading H-bond RGB color\n");
    if (fgets(str, sizeof str, fp) == NULL)  /* skip one line */
        fatal("error in reading comment line\n");
    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, "%lf %lf %lf", &width3[1], &width3[2], &width3[3]) != 3)
        fatal("error cylinder radius for bp-center line & bp 1 & 2\n");
    if (fgets(str, sizeof str, fp) == NULL)  /* skip one line */
        fatal("error in reading comment line\n");
    for (i = 0; i <= NATOMCOL; i++)
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &atom_col[i][1], &atom_col[i][2], &atom_col[i][3]) != 3)
            fatal("error reading atom RGB color\n");
    if (fgets(str, sizeof str, fp) == NULL)  /* skip one line */
        fatal("error in reading comment line\n");
    if (fgets(str, sizeof str, fp) == NULL)
        fatal("error reading label style\n");
    strcpy(label_style, str);

    close_file(fp);
}

/* write a record of round-ended cylinder (itype = 3) or flat-ended (5) or comments */
void r3d_rod(long itype, double *xyz1, double *xyz2, double rad, double *rgbv, FILE * fp)
{
    static char *format = "%9.3f";
    long i;

    if (itype == 3 || itype == 5)
        fprintf(fp, "%ld\n", itype);
    else
        fprintf(fp, "#5\n#");  /* comments and default to type 5 */
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz1[i]);
    fprintf(fp, format, rad);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz2[i]);
    fprintf(fp, format, rad);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, rgbv[i]);
    fprintf(fp, "\n");
}

/* draw dashed line for H-bonds, and lines connecting base-pair centers */
void r3d_dash(double *xyz1, double *xyz2, double hb_width, double *hb_col, FILE * fp)
{
    double distance, dnum, dxyz[4], m1[4], m2[4];
    long i, j, num;

    ddxyz(xyz1, xyz2, dxyz);
    distance = veclen(dxyz);
    vec_norm(dxyz);  /* normal vector */

    num = lround(hb_col[4] * distance);
    if (num % 2 == 0)
        num++;  /* odd "num" for ending on last point */
    if (num <= 1)  /* solid line */
        r3d_rod(5, xyz1, xyz2, hb_width, hb_col, fp);
    else {
        dnum = (double) num;
        for (i = 0; i < num; i += 2) {
            for (j = 1; j <= 3; j++) {
                m1[j] = xyz1[j] + dxyz[j] * distance * i / dnum;
                m2[j] = xyz1[j] + dxyz[j] * distance * (i + 1) / dnum;
            }
            r3d_rod(5, m1, m2, hb_width, hb_col, fp);
        }
    }
}

/* write a record of sphere */
void r3d_sphere(double *xyz1, double rad, double *rgbv, FILE * fp)
{
    static char *format = "%9.3f";
    long i;

    fprintf(fp, "2\n");
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz1[i]);
    fprintf(fp, format, rad);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, rgbv[i]);
    fprintf(fp, "\n");
}

void cpk_model(long num, long *idx, double **xyz, double ballrad, double **colrgb, FILE * fp)
{
    double vdw_radii[NELE];
    long i;

    if (ballrad <= 0.0)
        return;

    atom_info(3, NULL, NULL, vdw_radii);
    fprintf(fp, "###\n### The following section is for ball/CPK model\n");
    for (i = 1; i <= num; i++)
        r3d_sphere(xyz[i], vdw_radii[idx[i]] * ballrad, colrgb[i], fp);
}

/* write a record of triangle (itype = 1, default) or plane (itype = 6) */
void r3d_tripln(long itype, double *xyz1, double *xyz2, double *xyz3, double *rgbv, FILE * fp)
{
    static char *format = "%9.3f";
    long i;

    fprintf(fp, "%d\n", (itype != 6) ? 1 : 6);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz1[i]);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz2[i]);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz3[i]);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, rgbv[i]);
    fprintf(fp, "\n");
}

void r3d_block_edge(double *rgbv, long ioffset8, double **blkxyz, double w1, FILE * fp)
{
    static long blk_lkg[12][2] = {
        {1, 2},
        {2, 3},
        {3, 4},
        {1, 4},
        {1, 5},
        {2, 6},
        {3, 7},
        {4, 8},
        {5, 6},
        {6, 7},
        {7, 8},
        {5, 8}
    };
    long j;

    for (j = 0; j < 12; j++)
        r3d_rod(3, blkxyz[ioffset8 + blk_lkg[j][0]],
                blkxyz[ioffset8 + blk_lkg[j][1]], w1, rgbv, fp);
}

/* label the base in the center of six-membered ring */
void base_label(double **rxyz, char *label_style, double *rgbv, char *bname_num, FILE * fp)
{
    static char *format = "%9.3f";
    double cxyz[4];
    long i;

    avexyz(rxyz[1], rxyz[4], cxyz);  /* N1 + C4 */
    fprintf(fp, "10\n%s", label_style);  /* label_style has \n */
    fprintf(fp, "11\n");
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, cxyz[i]);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, rgbv[i]);
    fprintf(fp, "\n%s\n", bname_num);
}

void fill_base_ring(long num_residue, long num_ring, long **ring_atom, double **xyz,
                    long *ibase, char *bseq, double **base_col, char *label_style,
                    long label_ring, long *ResSeq, FILE * fp)
{
    char bname_num[BUF512];
    double **rxyz;
    long i, j;

    rxyz = dmatrix(1, 9, 1, 3);  /* base ring atoms */

    fprintf(fp, "###\n### The following section is for %ld filled base rings\n", num_ring);
    for (i = 1; i <= num_residue; i++) {
        if (ring_atom[i][10] <= 0)
            continue;
        for (j = 1; j <= ring_atom[i][10]; j++)
            cpxyz(xyz[ring_atom[i][j]], rxyz[j]);
        r3d_tripln(1, rxyz[3], rxyz[4], rxyz[6], base_col[ibase[i]], fp);
        r3d_tripln(1, rxyz[1], rxyz[2], rxyz[3], base_col[ibase[i]], fp);
        r3d_tripln(1, rxyz[4], rxyz[5], rxyz[6], base_col[ibase[i]], fp);
        r3d_tripln(1, rxyz[1], rxyz[3], rxyz[6], base_col[ibase[i]], fp);
        if (ring_atom[i][10] == 9) {
            r3d_tripln(1, rxyz[6], rxyz[7], rxyz[8], base_col[ibase[i]], fp);
            r3d_tripln(1, rxyz[1], rxyz[8], rxyz[9], base_col[ibase[i]], fp);
            r3d_tripln(1, rxyz[1], rxyz[6], rxyz[8], base_col[ibase[i]], fp);
        }
        if (label_ring) {
            sprintf(bname_num, "%c%ld", bseq[i], ResSeq[ring_atom[i][1]]);
            base_label(rxyz, label_style, base_col[ibase[i]], bname_num, fp);
        }
    }

    free_dmatrix(rxyz, 1, 9, 1, 3);
}

void process_alc(char *alcfile, char *imgfile, double scale_factor, long *opts)
{
    char **AtomName;
    double **blkxyz, **oxyz, **xyz;
    long i, j, k, nbond, nobj, npoint_blk, num, num_blk = 0, num_blk8;
    long ioffset8, ioffset_blk, nN = 0, nO = 0, nO_all = 0, nO_lkg = 0, nO_tmp;
    long Nidx[6], *blkibase, *ibase, **linkage;

    /* read in ALCHEMY file */
    get_alc_nums(alcfile, &num, &nbond);
    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 0, 2);
    read_alc(alcfile, &num, &nbond, AtomName, xyz, ibase, linkage);

    /* assume the first 8 atoms are 4 Ns followed by 4 Cs in
       a specific numbering scheme */
    for (i = 1; i <= num; i++)
        if (AtomName[i][0] == 'N')
            ++nN;
    if (!nN || nN % 4)
        fatal("wrong type of standard block file\n");
    else
        num_blk = nN / 4;

    nN = 0;
    for (i = 1; i <= num && nN < 5; i++)
        if (AtomName[i][0] == 'N')
            Nidx[++nN] = i;
    npoint_blk = (num_blk > 1) ? Nidx[5] - Nidx[1] : num;

    /* get the xyz coordinates of 8-point-per-block for each */
    num_blk8 = num_blk * 8;
    blkxyz = dmatrix(1, num_blk8, 1, 3);
    blkibase = lvector(1, num_blk);
    for (i = 1; i <= num_blk; i++) {
        ioffset8 = (i - 1) * 8;
        ioffset_blk = (i - 1) * npoint_blk;
        blkibase[i] = ibase[ioffset_blk + 1];
        if (!lval_in_range(blkibase[i], 0, 5)) {
            fprintf(stderr, "uncommon base. set to default color coding\n");
            blkibase[i] = NON_WC_IDX;
        }
        for (j = 1; j <= 8; j++)
            cpxyz(xyz[ioffset_blk + j], blkxyz[ioffset8 + j]);
    }
    nobj = num_blk * 6;  /* blocks */

    /* get the total number of O atoms:
       "O ", "OH", "OO", "OX", "OY", "OZ" */
    for (i = 1; i <= num; i++)
        if (AtomName[i][0] == 'O')
            nO_all++;
    oxyz = (nO_all) ? dmatrix(1, nO_all, 1, 3) : NULL;

    if (opts[3]) {  /* draw a line to link neighbor origins */
        nO_tmp = 0;
        for (i = 1; i <= num; i++)
            if (!strcmp(AtomName[i], "O ")) {
                ++nO_tmp;
                cpxyz(xyz[i], oxyz[nO_tmp]);
            }
        if (nO_tmp) {
            k = nO_tmp - 1;
            for (i = 1; i <= k; i++) {
                linkage[i][0] = 0;
                linkage[i][1] = i;
                linkage[i][2] = i + 1;
            }
            nO_lkg = k;
            nO = nO_tmp;
        }
    }
    if (opts[5]) {  /* draw helix-axis */
        nO_tmp = 0;
        for (i = 1; i <= num; i++)
            if (!strcmp(AtomName[i], "OH")) {
                ++nO_tmp;
                cpxyz(xyz[i], oxyz[nO + nO_tmp]);
            }
        if (nO_tmp) {
            if (nO_tmp % 2)
                fatal("wrong number of helix axis origins\n");
            k = nO_tmp / 2;
            for (i = 1; i <= k; i++) {
                j = nO_lkg + i;
                linkage[j][0] = -1;
                linkage[j][1] = nO + 2 * i - 1;
                linkage[j][2] = nO + 2 * i;
            }
            nO_lkg += k;
            nO += nO_tmp;
        }
    }
    if (opts[6]) {  /* draw global reference axis */
        nO_tmp = 0;
        for (i = 1; i <= num; i++)
            if (!strcmp(AtomName[i], "OO") || !strcmp(AtomName[i], "OX") ||
                !strcmp(AtomName[i], "OY") || !strcmp(AtomName[i], "OZ")) {
                ++nO_tmp;
                cpxyz(xyz[i], oxyz[nO + nO_tmp]);
            }
        if (nO_tmp) {
            k = nO_tmp - 1;
            for (i = 1; i <= k; i++) {
                j = nO_lkg + i;
                linkage[j][0] = -(i + 1);
                linkage[j][1] = nO + 1;
                linkage[j][2] = nO + i + 1;
            }
            nO_lkg += k;
            nO += nO_tmp;
        }
    }

    alc_3images(opts, nobj, num_blk, num_blk8, nO_lkg, nO, oxyz, blkxyz, blkibase,
                linkage, scale_factor, imgfile);

    if (oxyz != NULL)
        free_dmatrix(oxyz, 1, nO_all, 1, 3);
    free_alc(num, nbond, AtomName, xyz, ibase, 0, linkage);
    free_dmatrix(blkxyz, 1, num_blk8, 1, 3);
    free_lvector(blkibase, 1, num_blk);
}

/* From ALCHEMY format to three image formats: PS, XFIG, Raster3D */
void alc_3images(long *opts, long nobj, long num_blk, long num_blk8, long nO_lkg,
                 long nO, double **oxyz, double **blkxyz, long *blkibase, long **linkage,
                 double scale_factor, char *imgfile)
{
    char label_style[BUF512];
    double width3[4], hb_col[5], *pd, **atom_col, **base_col;
    static long faces[6][5] = {
        {1, 2, 3, 4, 1},  /* 0 front: minor groove */
        {5, 6, 7, 8, 5},  /* 1 back */
        {1, 4, 8, 5, 1},  /* 2 left */
        {2, 3, 7, 6, 2},  /* 3 right */
        {1, 2, 6, 5, 1},  /* 4 upper */
        {3, 4, 8, 7, 3}  /*  5 lower */
    };
    long default_size[2] = { PS_DFTSIZE, FIG_DFTSIZE };  /* PS & XFIG */
    long is_color, same_faces, updown, mgroove;
    long i, j, k, ioffset8, urxy[3];
    long *depth, *idx, **allobj;
    FILE *fp;

    is_color = opts[2];  /* decomposed for clarity: cf. alc2ps/alc2fig */
    same_faces = opts[7];
    updown = opts[9];
    mgroove = opts[10];

    fp = open_file(imgfile, "w");
    if (opts[0] == 2) {  /* Raster3D */
        atom_col = dmatrix(0, NATOMCOL, 1, 3);
        base_col = dmatrix(0, NBASECOL, 1, 3);
        raster3d_header(num_blk8, blkxyz, scale_factor, opts[8], opts[1], fp);
        /* only base_col is used here */
        get_r3dpars(base_col, hb_col, width3, atom_col, label_style);
        for (i = 1; i <= num_blk; i++) {
            ioffset8 = (i - 1) * 8;
            k = is_color ? blkibase[i] : NON_WC_IDX;  /* if color image */
            r3d_block_edge(same_faces ? base_col[NON_WC_IDX] : base_col[k],  /* if same-face */
                           ioffset8, blkxyz, width3[2], fp);  /* bp1 width */
            for (j = 0; j < 6; j++) {
                if (same_faces)
                    pd = (mgroove && !j) ? base_col[NON_WC_IDX] : base_col[k];
                else {
                    if ((!updown && !j) || (updown && j == 4))
                        pd = base_col[k];
                    else
                        pd = base_col[NBASECOL];
                }
                r3d_tripln(1, blkxyz[ioffset8 + faces[j][0]], blkxyz[ioffset8 + faces[j][1]], blkxyz[ioffset8 + faces[j][2]], pd, fp);  /* 0-1-2 */
                r3d_tripln(1, blkxyz[ioffset8 + faces[j][0]], blkxyz[ioffset8 + faces[j][2]], blkxyz[ioffset8 + faces[j][3]], pd, fp);  /* 0-2-3 */
            }
        }
        if (!is_color)
            cpxyz(base_col[NON_WC_IDX], hb_col);
        for (i = 1; i <= nO_lkg; i++) {
            if (!linkage[i][0])  /* bp-center connecting line */
                r3d_dash(oxyz[linkage[i][1]], oxyz[linkage[i][2]], width3[1], hb_col, fp);
            else  /* helix or reference frame */
                r3d_rod(3, oxyz[linkage[i][1]], oxyz[linkage[i][2]],
                        (linkage[i][0] == -1) ? width3[3] : width3[2], hb_col, fp);
        }

        free_dmatrix(atom_col, 0, NATOMCOL, 1, 3);
        free_dmatrix(base_col, 0, NBASECOL, 1, 3);
    } else {
        nobj += nO_lkg;  /* total number of objects */

        /* combine all objects (faces + origin lines) together */
        allobj = lmatrix(1, 3, 1, nobj);
        get_alc_objs(num_blk, blkxyz, nO, oxyz, nO_lkg, linkage, faces, allobj);

        idx = lvector(1, nobj);
        lsort(nobj, allobj[3], idx);

        adjust_xy(num_blk8, blkxyz, nO, oxyz, scale_factor, default_size[opts[0]], urxy);

        if (opts[0] == 1) {  /* for XFIG */
            depth = lvector(1, nobj);
            get_depth(nobj, allobj[3], depth);
            get_fig_xy(num_blk8, blkxyz, nO, oxyz, urxy, opts[1], fp);
            alc2fig(nobj, idx, depth, allobj, blkxyz, oxyz, blkibase, faces, opts, fp);
            free_lvector(depth, 1, nobj);
        } else {  /* PS */
            get_ps_xy(imgfile, urxy, opts[1], fp);
            alc2ps(nobj, idx, allobj, blkxyz, oxyz, blkibase, faces, opts, fp);
        }
        free_lmatrix(allobj, 1, 3, 1, nobj);
        free_lvector(idx, 1, nobj);
    }
    close_file(fp);
}

/* "allobj" has three row as follows
 *  (a) for blocks        [1] block index (1 -- num_blk)
 *                        [2] face index (0 -- 5)
 *                        [3] average z-coord x 1000 for each face
 *  (b) for O linkages    [1] 0: origin, -1: helix, -2, -3, -4: frame x, y, z
 *                        [2] 10000 x id1 + id2
 *                        [3] average z-coord x 1000 for the line */
void get_alc_objs(long num_blk, double **blkxyz, long nO, double **oxyz, long nO_lkg,
                  long **linkage, long faces[][5], long **allobj)
{
    long i, ioffset8, ip, j, k;

    for (i = 1; i <= num_blk; i++) {
        ioffset8 = (i - 1) * 8;  /* ith block atom index offset */
        k = (i - 1) * 6 + 1;  /* ith block face offset */
        for (j = 0; j < 6; j++) {
            ip = k + j;
            allobj[1][ip] = i;
            allobj[2][ip] = j;
            /* only two diagonal atoms are needed for z-coordinates */
            allobj[3][ip] = lround(500.0 * (blkxyz[ioffset8 + faces[j][0]][3] +
                                            blkxyz[ioffset8 + faces[j][2]][3]));
        }
    }

    if (nO) {
        k = num_blk * 6;  /* current object offset */
        for (i = 1; i <= nO_lkg; i++) {
            j = k + i;
            allobj[1][j] = linkage[i][0];
            allobj[2][j] = linkage[i][1] * 10000 + linkage[i][2];
            allobj[3][j] = lround(500.0 * (oxyz[linkage[i][1]][3] + oxyz[linkage[i][2]][3]));
        }
    }
}

/* read in rendering parameters for XFIG file */
void get_fig_pars(double *dot_sep, long *dlcol, long *dwidth, long *bp1width,
                  long *bp2width, long **bc_idx, double *msat, double *Msat,
                  long *o_sides, long *line_width, long *join_style, long *cap_style, long *mfcol)
{
    char BDIR[BUF512], str[BUF512];
    char *fig_image_par = "fig_image.par";
    long i;
    FILE *fp;

    get_BDIR(BDIR, fig_image_par);
    strcat(BDIR, fig_image_par);
    fprintf(stderr, " ...... reading file: %s ...... \n", fig_image_par);

    fp = open_file(BDIR, "r");
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%lf", dot_sep) != 1)
        fatal("error in reading dot separation\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", dlcol) != 1)
        fatal("error in reading dot line color\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", dwidth) != 1)
        fatal("error in reading dot line width\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", bp1width) != 1)
        fatal("error in reading bp1 width\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", bp2width) != 1)
        fatal("error in reading bp2 width\n");

    for (i = 0; i <= 6; i++)
        if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", &bc_idx[1][i]) != 1)
            fatal("error in reading color-code\n");
    for (i = 0; i <= 6; i++)  /* black & white option */
        bc_idx[0][i] = 0;

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%lf", msat) != 1)
        fatal("error in reading minor groove color saturation\n");
    *msat = fabs(*msat);
    if (*msat > 1.0)
        *msat = 0.9;
    *msat = 1.0 - *msat;
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%lf", Msat) != 1)
        fatal("error in reading major groove color saturation\n");
    *Msat = fabs(*Msat);
    if (*Msat > 1.0)
        *Msat = 0.1;
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", o_sides) != 1)
        fatal("error in reading other-side color code\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", line_width) != 1)
        fatal("error in reading line width\n");
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", join_style) != 1)
        fatal("error in reading join style\n");
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", cap_style) != 1)
        fatal("error in reading cap style\n");
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", mfcol) != 1)
        fatal("error in reading minor groove filling color\n");

    close_file(fp);
}

/* align the structure so that its global reference frame is in that
 * defined in 'morg' and 'mst' (column-wise for x-, y- and z-axes) */
void frame_xyz(long side_view, double *morg, double **mst, long num, double **xyz)
{
    double **tmpxyz;
    long i, j;

    tmpxyz = dmatrix(1, num, 1, 3);

    for (i = 1; i <= num; i++)
        for (j = 1; j <= 3; j++)
            tmpxyz[i][j] = dot(xyz[i], mst[j]) + morg[j];

    copy_dmatrix(tmpxyz, num, 3, xyz);

    if (side_view)
        get_side_view(1, num, xyz);

    free_dmatrix(tmpxyz, 1, num, 1, 3);
}

void change_xyz(long side_view, double *morg, double **mst, long num, double **xyz)
{
    double **tmpxyz;

    tmpxyz = dmatrix(1, num, 1, 3);

    move_position(xyz, num, 3, morg);
    multi_matrix(xyz, num, 3, mst, 3, 3, tmpxyz);

    copy_dmatrix(tmpxyz, num, 3, xyz);

    if (side_view)
        get_side_view(1, num, xyz);

    free_dmatrix(tmpxyz, 1, num, 1, 3);
}

/* adjust orientation: xyz*rotz(-90)*rotx(90) ==> [-y z -x] */
void get_side_view(long ib, long ie, double **xyz)
{
    long i;
    double temp;

    for (i = ib; i <= ie; i++) {
        temp = xyz[i][1];
        xyz[i][1] = -xyz[i][2];
        xyz[i][2] = xyz[i][3];
        xyz[i][3] = -temp;
    }
}

/* get C1* and RN9/YN1 atom pair index for finding helix */
void get_CNidx(long ds, long num_bp, long **chi, long **idx, long *nvec, long *C1b, long *C1e)
{
    long i, ia, ib, ioffset, j, joffset, k;

    *nvec = 0;
    for (i = 1; i <= ds; i++)
        for (j = 1; j <= num_bp - 1; j++) {
            ioffset = (j - 1) * 4;
            joffset = j * 4;
            for (k = 2; k <= 3; k++) {
                ia = chi[i][ioffset + k];
                ib = chi[i][joffset + k];
                if (ia && ib) {
                    idx[++*nvec][1] = ia;
                    idx[*nvec][2] = ib;
                }
            }
        }

    ioffset = (num_bp - 1) * 4;
    for (i = 1; i <= ds; i++) {
        C1b[i] = chi[i][2];
        C1e[i] = chi[i][ioffset + 2];
    }
}

/* attach reference frame axes */
void add_3axes(long *num, char **AtomName, long *ibase, double **xyz, long *nbond,
               long **linkage, long side_view, double axis_len)
{
    static char *ref_symbol[4] = { "OO", "OX", "OY", "OZ" };
    static double ref_xyz[4][3] = {
        {0.0, 0.0, 0.0},  /* origin */
        {1.0, 0.0, 0.0},  /* x-axis */
        {0.0, 1.0, 0.0},  /* y-axis */
        {0.0, 0.0, 1.0}  /* z-axis */
    };
    long i, j, k;

    for (i = 1; i <= 4; i++) {
        k = *num + i;
        strcpy(AtomName[k], ref_symbol[i - 1]);
        for (j = 1; j <= 3; j++)
            xyz[k][j] = axis_len * ref_xyz[i - 1][j - 1];
        ibase[k] = NON_WC_IDX;  /* non-common base atoms */
    }
    if (side_view)
        get_side_view(*num + 1, *num + 4, xyz);
    for (i = 1; i <= 3; i++) {
        k = *nbond + i;
        linkage[k][1] = *num + 1;
        linkage[k][2] = *num + i + 1;
    }
    *num += 4;
    *nbond += 3;
}

char *get_sequence(char *Wbase, long *num_bp)
{
    char *bseq;
    char str[BUF512];
    long ich = 2, k, nbp;

    fprintf(stderr, "\nInput your base sequence with only %s:\n", Wbase);
    fprintf(stderr, "1. From a data file (complete sequence)\n");
    fprintf(stderr, "2. From keyboard (enter only the repeating sequence)\n");
    fprintf(stderr, "Your choice (1 or 2, Dft: 2): ");
    fflush(stderr);
    if (fgets(str, sizeof str, stdin) != NULL) {
        k = sscanf(str, "%ld", &ich);
        if (!k || k == EOF)
            ich = 2;
    } else
        fatal("error in reading your choice\n");
    fprintf(stderr, "\n");

    if (ich == 1) {
        fprintf(stderr, "Name of your base sequence file: ");
        fflush(stderr);
        if (fgets(str, sizeof str, stdin) == NULL)
            fatal("error in reading your sequence file name\n");
        trim(str);
        bseq = read_sequence(str, Wbase, &nbp);
    } else
        bseq = read_repeat("A", 0, Wbase, &nbp);

    *num_bp = nbp;

    return bseq;
}

char **single2double(long nbp, char *bseq, char *Wbase, char *Cbase)
{
    char **bp_seq;
    long i, j;

    bp_seq = cmatrix(1, nbp, 1, 2);
    for (i = 1; i <= nbp; i++) {
        bp_seq[i][1] = bseq[i - 1];
        j = strchr(Wbase, bseq[i - 1]) - Wbase;
        bp_seq[i][2] = *(Cbase + j);
    }
    free_cvector(bseq, 0, nbp);

    return bp_seq;
}

long is_valid_base(char c, char *valid_bases)
{
    if (isspace((int) c))  /* skip white space */
        return FALSE;

    if (strchr(valid_bases, c) == NULL) {
        fprintf(stderr, "skip %c: acceptable bases [%s]\n", c, valid_bases);
        return FALSE;
    }

    return TRUE;
}

long repeat_num(void)
{
    char *p0;
    long num = 10, k;

    fprintf(stderr, "Number of repeats (Dft: 10): ");
    fflush(stderr);
    if ((p0 = my_getline(stdin)) != NULL) {
        k = sscanf(p0, "%ld", &num);
        if (!k || k == EOF || !num)
            num = 10;
    } else
        fatal("error in reading number of repeats\n");
    free(p0);

    if (num < 0) {
        fprintf(stderr, "number of repeat %ld < 0. reset to positive\n", num);
        num = -num;
    }

    return num;
}

/* read in sequence information from file 'seqfile', return number of
 * valid bases [must be in UPPER case], and array of the bases */
char *read_sequence(char *seqfile, char *valid_bases, long *nbp)
{
    char *p0, *line, *pb, *bseq;
    long nb = 0, maxline = BUF512;
    FILE *fp;

    pb = (char *) malloc(maxline * sizeof(char));  /* initial size */
    if (pb == NULL)
        fatal("malloc failure in reading sequence\n");

    fp = open_file(seqfile, "r");
    while ((p0 = my_getline(fp)) != NULL) {
        line = trim(p0);  /* keep the original value of p0 */
        if (line[0] == '#' || line[0] == '\0') {
            free(p0);
            continue;  /* skip empty and commented lines */
        }
        upperstr(line);  /* changed to upper case */
        while (*line) {
            if (is_valid_base(*line, valid_bases)) {
                if (nb >= maxline - 1)
                    pb = enlarge_cline(&maxline, pb);
                pb[nb++] = *line;
            }
            line++;
        }
        free(p0);
    }
    pb[nb] = '\0';

    close_file(fp);
    bseq = my_strdup(pb);
    free(pb);

    if (!nb)
        fatal("sequence file <%s> contains no valid bases\n", seqfile);

    *nbp = nb;

    return bseq;
}

/* read in repeated sequence, return number of valid bases [must be in
 * UPPER case], and array of the bases */
char *read_repeat(char *crepeat, long fixed, char *valid_bases, long *nbp)
{
    char *p0, *line, *pb, *bseq;
    long i, k, num, num_total, nb = 0, maxline = BUF512;

    pb = (char *) malloc(maxline * sizeof(char));
    if (pb == NULL)
        fatal("malloc failure in reading sequence\n");

    if (fixed) {
        fprintf(stderr, "Repeating unit: %s\n", crepeat);
        fflush(stderr);
        strcpy(pb, crepeat);
        nb = strlen(crepeat);
    } else {
        fprintf(stderr, "Repeating unit (Dft: %s): ", crepeat);
        fflush(stderr);
        if ((p0 = my_getline(stdin)) == NULL)
            fatal("error in reading your repeating unit\n");
        line = trim(p0);
        upperstr(line);

        while (*line) {
            if (is_valid_base(*line, valid_bases)) {
                if (nb >= maxline - 1)
                    pb = enlarge_cline(&maxline, pb);
                pb[nb++] = *line;
            }
            line++;
        }
        free(p0);

        if (!nb) {  /* using default */
            strcpy(pb, crepeat);
            nb = strlen(crepeat);
        } else
            pb[nb] = '\0';

        fprintf(stderr, "Repeating unit: %s\n", pb);
    }

    num = repeat_num();

    num_total = nb * num;
    bseq = cvector(0, num_total);

    for (i = 1; i <= num; i++) {
        k = (i - 1) * nb;
        strcpy(bseq + k, pb);
    }

    free(pb);

    *nbp = num_total;

    return bseq;
}

/* combine strand II with strand I for a parallel duplex */
void combine_pstnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                    char **tAtomName, char **tResName, char *tChainID, long *tResSeq,
                    double **txyz, char **tAtomName2, double **txyz2)
{
    long i, ik, j;

    for (i = 1; i <= num_bp; i++) {
        for (j = s2idx[i][1] + 1; j <= s2idx[i][1] + s2idx[i][2]; j++) {
            ik = *tnum + j - s2idx[i][1];
            strcpy(tAtomName[ik], tAtomName2[j]);
            sprintf(tResName[ik], "  %c", bp_seq[i][2]);
            tChainID[ik] = Gvars.REBUILD_CHAIN_IDS[1];
            tResSeq[ik] = num_bp + i;
            cpxyz(txyz2[j], txyz[ik]);
        }
        *tnum += s2idx[i][2];
    }
}

/* reverse strand II to be in 5'-->3' direction & combined with I */
void reverse_stnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                   char **tAtomName, char **tResName, char *tChainID, long *tResSeq,
                   double **txyz, char **tAtomName2, double **txyz2, long basep)
{
    long i, ik, j;

    for (i = num_bp; i >= 1; i--) {
        for (j = s2idx[i][1] + 1; j <= s2idx[i][1] + s2idx[i][2]; j++) {
            ik = *tnum + j - s2idx[i][1];
            strcpy(tAtomName[ik], tAtomName2[j]);
            sprintf(tResName[ik], "  %c", bp_seq[i][2]);
            tChainID[ik] = Gvars.REBUILD_CHAIN_IDS[1];
            tResSeq[ik] = 2 * num_bp - i + 1;
            if (basep && !strcmp(tAtomName[ik], " P  "))
                ++tResSeq[ik];  /* re-position P on strand II */
            cpxyz(txyz2[j], txyz[ik]);
        }
        *tnum += s2idx[i][2];
    }
}

void pair_checking(long ip, long ds, long num_residue, char *pdbfile, long *num_bp,
                   long **pair_num)
{
    long i, j;

    if (!ip) {  /* ideal base-pairing */
        if (ds == 1) {
            if (*num_bp > num_residue) {
                *num_bp = num_residue;
                fprintf(stderr, "processing all %ld residues\n\n", *num_bp);
            }
        } else if (num_residue % 2) {  /* ds = 2 */
            fprintf(stderr, "%s has odd number (%ld) of residues\n", pdbfile, num_residue);
            fatal("Please specify base-pairing residue numbers\n");
        } else {
            i = num_residue / 2;
            if (*num_bp < i)
                fatal("Please specify base-pairing residue numbers\n");
            if (*num_bp > i) {
                *num_bp = i;
                fprintf(stderr, "processing all %ld base-pairs\n\n", i);
            }
        }

        for (i = 1; i <= *num_bp; i++) {
            pair_num[1][i] = i;
            if (ds == 2)
                pair_num[2][i] = 2 * *num_bp - i + 1;
        }
    } else {
        if (ds * *num_bp > num_residue)
            fprintf(stderr, "some residue has more than one pair\n");
        for (i = 1; i <= ds; i++)
            for (j = 1; j <= *num_bp; j++)
                if (pair_num[i][j] > num_residue) {
                    fprintf(stderr, "residue index %ld too big (> %ld)\n",
                            pair_num[i][j], num_residue);
                    fatal("please check your input file\n");
                }
    }
}

void double_print_msg(char *msg, FILE * fp)
{
    fprintf(stderr, "%s\n", msg);
    fprintf(fp, "%s\n", msg);
}

/* check for proper O3'-P linkage, strand I 5'--> 3' direction, parallel duplex */
void drct_checking(long ds, long num_bp, long **pair_num, long **seidx, char **AtomName,
                   double **xyz, long *parallel, long *bbexist, long **o3p_brk, FILE * fp)
{
    char msg[BUF512];
    long i, ib, ie, ioffset, j, jm1, jp1, rnum, n = 0;
    long direction[7], **O3, **P;

    O3 = lmatrix(1, 2, 1, num_bp);
    P = lmatrix(1, 2, 1, num_bp);

    for (i = 1; i <= ds; i++)
        for (j = 1; j <= num_bp; j++) {
            rnum = pair_num[i][j];
            ib = seidx[rnum][1];
            ie = seidx[rnum][2];
            P[i][j] = find_1st_atom(" P  ", AtomName, ib, ie, "");
            O3[i][j] = find_1st_atom(" O3'", AtomName, ib, ie, "");
            if (P[i][j] && O3[i][j]) {
                ++*bbexist;
                if (within_limits(xyz[O3[i][j]], xyz[P[i][j]], 0.8, O3P_UPPER))
                    n++;
            }
        }

    if (n) {
        print_sep(fp, '*', 76);
        sprintf(msg, "WARNING: %ld out of %ld bases have O3'[i] wrongly connected"
                " to P[i]", n, ds * num_bp);
        double_print_msg(msg, fp);
    }

    init_lvector(direction, 1, 6, 0);
    for (i = 1; i <= ds; i++) {
        ioffset = (i - 1) * 3;
        for (j = 1; j < num_bp; j++) {
            jp1 = j + 1;
            if (O3[i][j] && P[i][jp1] &&  /* j to j + 1 */
                within_limits(xyz[P[i][jp1]], xyz[O3[i][j]], 0.8, O3P_UPPER)) {
                ++direction[ioffset + 1];
                continue;
            }
            if (O3[i][jp1] && P[i][j] &&  /* j + 1 to j */
                within_limits(xyz[P[i][j]], xyz[O3[i][jp1]], 0.8, O3P_UPPER)) {
                ++direction[ioffset + 2];
                continue;
            }
            ++direction[ioffset + 3];
            if (*bbexist)
                o3p_brk[i][j] = 9;
        }
    }

    if ((direction[1] && direction[2]) || (direction[4] && direction[5])) {
        sprintf(msg, "This structure contains intra-chain direction reverse");
        double_print_msg(msg, fp);
    } else {
        if (*bbexist && (direction[3] || direction[6])) {
            sprintf(msg, "This structure has broken O3'[i] to P[i+1] linkages");
            double_print_msg(msg, fp);
        }
        if (direction[2] && !direction[1]) {
            sprintf(msg, "WARNING: Strand I in 3'-->5' direction!!!\n");
            double_print_msg(msg, fp);
        }
        if (direction[4] && !direction[5]) {
            fprintf(stderr, "This is a parallel duplex\n");
            *parallel = 1;
        }
    }
    print_sep(fp, '*', 76);

    if (ds == 1) {
        long ka, kb, ksum;

        for (j = 1; j <= num_bp; j++) {
            ksum = 0;

            jm1 = j - 1;
            if (jm1 < 1)
                ksum++;
            else {
                ka = O3[1][jm1];
                kb = P[1][j];
                if (ka && kb && within_limits(xyz[ka], xyz[kb], 0.8, O3P_UPPER))
                    continue;
                ksum++;
            }

            jp1 = j + 1;
            if (jp1 > num_bp)
                ksum++;
            else {
                ka = O3[1][j];
                kb = P[1][jp1];
                if (ka && kb && within_limits(xyz[ka], xyz[kb], 0.8, O3P_UPPER))
                    continue;
                ksum++;
            }

            if (*bbexist && ksum == 2)
                o3p_brk[1][j] = 1;
        }
    }

    free_lmatrix(O3, 1, ds, 1, num_bp);
    free_lmatrix(P, 1, ds, 1, num_bp);
}

void residue_idstr(char chain_id, long res_seq, char *rname, char *idmsg)
{
    long i;

    sprintf(idmsg, "%c:%ld:%s", chain_id, res_seq, rname);
    for (i = 0; i < (long) strlen(idmsg); i++)
        if (idmsg[i] == ' ')
            idmsg[i] = '_';
}

/* get base name: all information to uniquely identify a base residue */
void base_str(char chain_id, long res_seq, char *misc, char *rname, char bcode,
              long stnd, char *idmsg)
{
    char b1, snum[10], rname_cp[10], modelNum[10], iCode = misc[2];
    long i, nlen = 4;

    b1 = (bcode == '\0') ? ' ' : bcode;

    strncpy(modelNum, misc + 30, nlen);
    modelNum[nlen] = '\0';
    for (i = 0; i < nlen; i++)
        if (modelNum[i] == '0')
            modelNum[i] = '.';
    sprintf(snum, "%4ld", res_seq);
    for (i = 0; i < 4; i++)
        if (snum[i] == ' ')
            snum[i] = '.';
    if (chain_id == ' ')
        chain_id = '-';
    if (iCode == ' ')
        iCode = '_';
    strcpy(rname_cp, rname);
    for (i = 0; i < 3; i++)
        if (rname_cp[i] == ' ')
            rname_cp[i] = '.';
    if (stnd == 1)  /* strand I */
        sprintf(idmsg, "%4s>%c:%4s%c:[%s]%c", modelNum, chain_id, snum, iCode, rname_cp, b1);
    else  /* strand II */
        sprintf(idmsg, "%c[%s]:%4s%c:%c<%4s", b1, rname_cp, snum, iCode, chain_id, modelNum);
}

/* write out a bond list for checking and easy modification */
void write_lkglist(long nbond, long **linkage, char **AtomName, char **ResName,
                   char *ChainID, long *ResSeq, char **Miscs)
{
    char b1[BUF512], b2[BUF512];
    long i, ia, ib;
    FILE *fp;

    fp = open_file(LKG_FILE, "w");
    for (i = 1; i <= nbond; i++) {
        ia = linkage[i][1];
        ib = linkage[i][2];
        base_str(ChainID[ia], ResSeq[ia], Miscs[ia], ResName[ia], '\0', 1, b1);
        base_str(ChainID[ib], ResSeq[ib], Miscs[ib], ResName[ib], '\0', 1, b2);
        fprintf(fp, "%5ld %5ld  # %5ld %s {%s} with %s {%s}\n", ia, ib, i,
                AtomName[ia], b1, AtomName[ib], b2);
    }
    close_file(fp);
}

/* get all possible H-bonds with upper/lower layer information */
void hbond_info(long **pair_num, char *bseq, long **seidx, long *idx, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                long *RY, long *num_hbond, long **hb_linkage)
{
    long i, j, k, ki, kj, num_bp = 2;
    long **num_list;
    FILE *fph;

    num_list = lmatrix(1, BUF512, 0, 3);
    fph = open_file(HB_FILE, "w");
    for (k = 1; k <= num_bp; k++) {
        for (i = 1; i < pair_num[k][0]; i++) {
            ki = pair_num[k][i];
            if (RY[ki] < 0)
                continue;
            for (j = i + 1; j <= pair_num[k][0]; j++) {  /* pairwise */
                kj = pair_num[k][j];
                if (RY[kj] < 0)
                    continue;
                hbond_list(ki, kj, AtomName, ResName, ChainID, ResSeq, xyz, Miscs,
                           bseq, seidx, idx, hb_linkage, &Gvars.misc_pars, num_list,
                           num_hbond, k, fph);
            }
        }
    }
    close_file(fph);
    free_lmatrix(num_list, 1, BUF512, 0, 3);
}

/* get H-bonds betwen base-pairs for a general PDB data file */
void hbond_pdb(long num, long num_residue, char *bseq, long **seidx, long *idx, long *RY,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, long *num_hbond, long **hb_linkage, long pwise)
{
    double rtn_val[RTNNUM];
    double **orien, **org, **NC1xyz, **o3_p;
    long bpid, i, j, **num_list, **ring_atom;
    FILE *fph;

    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 3);
    NC1xyz = dmatrix(1, num_residue, 1, 7);  /* RN9/YN1 & C1' atomic coordinates */
    o3_p = dmatrix(1, num_residue, 1, 8);  /* O3'/P atomic coordinates */

    base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID, ResSeq,
              Miscs, xyz, orien, org, NC1xyz, o3_p);

    /* 1-9 ring atom index, 10 # of ring atoms, 11-19 first level */
    ring_atom = lmatrix(1, num_residue, 1, 19);
    ring_oidx(num, num_residue, RY, seidx, AtomName, xyz, idx, ring_atom);

    num_list = lmatrix(1, BUF512, 0, 3);

    fph = open_file(HB_FILE, "w");
    for (i = 1; i < num_residue; i++) {
        if (!pwise && RY[i] < 0)
            continue;
        for (j = i + 1; j <= num_residue; j++) {
            if (!pwise) {  /* check for base pairing */
                if (RY[j] < 0)  /* only considering hydration for the 1st layer */
                    continue;
                check_pair(i, j, bseq, seidx, xyz, NC1xyz, orien, org, idx,
                           AtomName, &Gvars.misc_pars, rtn_val, &bpid, ring_atom, 0);
                if (!bpid)
                    continue;  /* not a base-pair */
            }
            hbond_list(i, j, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, bseq, seidx,
                       idx, hb_linkage, &Gvars.misc_pars, num_list, num_hbond, 1, fph);
        }
    }

    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 3);
    free_dmatrix(NC1xyz, 1, num_residue, 1, 7);
    free_dmatrix(o3_p, 1, num_residue, 1, 8);
    free_lmatrix(ring_atom, 1, num_residue, 1, 19);
    free_lmatrix(num_list, 1, BUF512, 0, 3);
    close_file(fph);
}

/* write out h-bond listing information: <num_hbond> already a pointer */
void hbond_list(long i, long j, char **AtomName, char **ResName, char *ChainID,
                long *ResSeq, double **xyz, char **Miscs, char *bseq, long **seidx,
                long *idx, long **hb_linkage, miscPars * misc_pars, long **num_list,
                long *num_hbond, long ilayer, FILE * fph)
{
    char b1[BUF512], b2[BUF512];
    long ir, jr, k, num_hb;

    ir = seidx[i][1];
    jr = seidx[j][1];
    base_str(ChainID[ir], ResSeq[ir], Miscs[ir], ResName[ir], bseq[i], 1, b1);
    base_str(ChainID[jr], ResSeq[jr], Miscs[jr], ResName[jr], bseq[j], 1, b2);
    hb_numlist(i, j, bseq[i], bseq[j], seidx, idx, AtomName, xyz, misc_pars, &num_hb, num_list);
    for (k = 1; k <= num_hb; k++) {
        hb_linkage[++*num_hbond][1] = ilayer;
        hb_linkage[*num_hbond][2] = num_list[k][1];
        hb_linkage[*num_hbond][3] = num_list[k][2];
        fprintf(fph, "%5ld %5ld  # %5ld %4.2f%c", num_list[k][1], num_list[k][2],
                *num_hbond, num_list[k][3] / MFACTOR, num_list[k][0] ? '*' : ' ');
        fprintf(fph, "%s {%s} with %s {%s}\n", AtomName[num_list[k][1]], b1,
                AtomName[num_list[k][2]], b2);
    }
}

/* determine atom number list for each base-pair from string hb_info */
void hb_numlist(long i, long j, char basei, char basej, long **seidx, long *idx,
                char **AtomName, double **xyz, miscPars * misc_pars, long *num_hb, long **num_list)
{
    char atom1[10], atom2[10], hb_info[BUF512], *pchar;
    long k, num;

    get_hbond_ij(i, j, basei, basej, misc_pars, seidx, idx, AtomName, xyz, hb_info);

    if (sscanf(hb_info, "[%ld]", &num) != 1)
        num = 0;
    if (num) {
        pchar = strchr(hb_info, ' ');
        for (k = 1; k <= num; k++) {
            strncpy(atom1, pchar + 1, 4);
            atom1[4] = '\0';
            strncpy(atom2, pchar + 6, 4);
            atom2[4] = '\0';
            num_list[k][1] = find_1st_atom(atom1, AtomName, seidx[i][1], seidx[i][2], "");
            num_list[k][2] = find_1st_atom(atom2, AtomName, seidx[j][1], seidx[j][2], "");

            num_list[k][0] = 0;
            if (*(pchar + 5) == '*')
                num_list[k][0] = 1;  /* not a canonical D-A atom pair */
            strncpy(atom1, pchar + 11, 4);  /* HB distance */
            atom1[4] = '\0';
            num_list[k][3] = lround(cvt2double(atom1) * MFACTOR);
            pchar += 15;  /* point to next H-bond */
        }
    }
    *num_hb = num;
}

/* get detailed H-bonding information for each pair */
void hb_information(long num_bp, long **pair_num, char **bp_seq, long **seidx, long *idx,
                    char **AtomName, double **xyz, long *WC_info, FILE * fp)
{
    char hb_info[BUF512];
    long i, j, k;

    print_sep(fp, '*', 76);
    fprintf(fp, "Detailed H-bond information: atom-name pair and length [%s]\n",
            Gvars.misc_pars.hb_atoms);

    for (k = 1; k <= num_bp; k++) {
        i = pair_num[1][k];
        j = pair_num[2][k];
        get_hbond_ij(i, j, bp_seq[1][k], bp_seq[2][k], &Gvars.misc_pars, seidx, idx,
                     AtomName, xyz, hb_info);
        fprintf(fp, "%4ld %c-%c%c%c-%c  %s\n", k, bp_seq[1][k],
                (WC_info[k] == 2) ? '-' : '*', WC_info[k] ? '-' : '*', bp_seq[0][k],
                bp_seq[2][k], hb_info);
    }
}

/* check if two atom pair could form a H-bond */
long good_hbatoms(miscPars * misc_pars, char *atom1, char *atom2, long idx1, long idx2)
{
    static char *PO[] = { " O1P", " O2P", " O3'", " O4'", " O5'", " N7 " };
    static long numPO = sizeof PO / sizeof PO[0] - 1;
    long natom = misc_pars->hb_idx[0];

    if (num_strmatch(atom1, PO, 0, numPO) && num_strmatch(atom2, PO, 0, numPO))
        return FALSE;  /* no H-bond within PO4 group, O4' & N7 */

    if ((idx1 == 2 || idx1 == 4 || idx2 == 2 || idx2 == 4) &&  /* at least one O/N atom */
        (lval_in_set(idx1, 1, natom, misc_pars->hb_idx) &&
         lval_in_set(idx2, 1, natom, misc_pars->hb_idx)))
        return TRUE;
    else
        return FALSE;
}

/* read in bond linkage information from "lkgfile": r3d_atom & stack2img generate
 * a file called "bonds_lkg.dat", which can be modified to fit any purpose */
void read_lkginfo(char *lkgfile, long num, long *nbond, long **linkage)
{
    char str[BUF512];
    long ia, ib;
    FILE *fp;

    fp = open_file(lkgfile, "r");
    while (fgets(str, sizeof str, fp) != NULL)
        if (sscanf(str, "%ld %ld", &ia, &ib) == 2) {
            if (!lval_in_range(ia, 1, num) || !lval_in_range(ib, 1, num)) {
                fprintf(stderr, "%s", str);
                fprintf(stderr, "         has atom serial number out of range\n");
                continue;
            }
            linkage[++*nbond][1] = ia;
            linkage[*nbond][2] = ib;
        }
    close_file(fp);
}

/* read in H-bonding information from "hbfile": r3d_atom & stack2img generate
 * a file called "hbonds_info.dat", which can be modified to fit any purpose */
void read_hbinfo(char *hbfile, long num, long *num_hbond, long **hb_linkage)
{
    char str[BUF512];
    long ia, ib;
    FILE *fp;

    fp = open_file(hbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL)
        if (sscanf(str, "%ld %ld", &ia, &ib) == 2) {
            if (!lval_in_range(ia, 1, num) || !lval_in_range(ib, 1, num)) {
                fprintf(stderr, "[%s] has atom serial number out of range\n", str);
                continue;
            }
            hb_linkage[++*num_hbond][1] = 1;
            hb_linkage[*num_hbond][2] = ia;
            hb_linkage[*num_hbond][3] = ib;
        }
    close_file(fp);
}

void update_hb_idx(long idx, double *dtmp, long *ddidx, double *hb_dist, long cur_idx)
{
    dtmp[idx] = hb_dist[cur_idx];
    ddidx[idx] = cur_idx;
}

/* get the H-bonding atom-pairs */
void hb_atompair(long num_hbonds, char **hb_atom1, char **hb_atom2, double *hb_dist,
                 long *lkg_type, miscPars * misc_pars)
{
    double dtmp[3];
    long k, m = 0, n, num_iter = 1;
    long ddidx[3], *matched_idx, **idx2;

    if (!num_hbonds)
        return;

    matched_idx = lvector(1, num_hbonds);

    while (1) {
        if (matched_idx[num_iter]) {
            num_iter++;
            continue;
        }

        for (k = 1; k <= 2; k++)
            update_hb_idx(k, dtmp, ddidx, hb_dist, num_iter);

        for (n = 1; n <= num_hbonds; n++) {
            if (n == num_iter || matched_idx[n])
                continue;
            if (!strcmp(hb_atom1[n], hb_atom1[num_iter]) && hb_dist[n] < dtmp[1])
                update_hb_idx(1, dtmp, ddidx, hb_dist, n);
            if (!strcmp(hb_atom2[n], hb_atom2[num_iter]) && hb_dist[n] < dtmp[2])
                update_hb_idx(2, dtmp, ddidx, hb_dist, n);
        }

        if (ddidx[1] == ddidx[2]) {  /* best mutual match */
            k = ddidx[1];
            hb_dist[k] = -hb_dist[k];  /* make  */
            num_iter = 1;  /* reset the iterator */

            for (n = 1; n <= num_hbonds; n++) {
                if (matched_idx[n])
                    continue;
                if (!strcmp(hb_atom1[n], hb_atom1[k]) ||  /* not && */
                    !strcmp(hb_atom2[n], hb_atom2[k])) {
                    matched_idx[n] = 1;
                    m++;
                }
            }

            if (m >= num_hbonds)
                break;

        } else
            num_iter++;
    }

    /* === further processing by adding more possible H-bonds === */
    idx2 = lmatrix(1, num_hbonds, 1, 2);

    for (k = 1; k <= num_hbonds; k++) {
        if (hb_dist[k] > 0.0)
            continue;
        idx2[k][1] = 9;
        idx2[k][2] = 9;
        for (m = 1; m <= num_hbonds; m++) {
            if (m == k || hb_dist[m] < 0.0)
                continue;
            if (!strcmp(hb_atom1[m], hb_atom1[k]))
                idx2[m][1] = 1;
            if (!strcmp(hb_atom2[m], hb_atom2[k]))
                idx2[m][2] = 1;
        }
    }

    /* Note that at this point, there are only 3 possibilities:
       9 + 9 = 18 means already located H-bonds;
       1 + 1 = 2 means the two atoms in separate H-bonds;
       1 + 0 = 0 + 1 = 1 means only one of atoms in H-bonds */
    for (k = 1; k <= num_hbonds; k++) {
        m = idx2[k][1] + idx2[k][2];
        lkg_type[k] = m;  /* for later on referencing */
        if (m != 18 && dval_in_range(hb_dist[k], misc_pars->hb_lower, misc_pars->hb_dist2))
            hb_dist[k] = -hb_dist[k];
    }

    /* The following code for adding additional H-bonds is too
     * complicated, ad hoc, but seems adding not much useful information
     * as from the above simple approach -- feb28, 2008 */

/* ============================================================================ */
#if 0
    for (k = 1; k <= num_hbonds; k++) {
        m = idx2[k][1] + idx2[k][2];
        lkg_type[k] = m;  /* for later on referencing */
        if (m == 1 && dval_in_range(hb_dist[k], misc_pars->hb_lower, misc_pars->hb_dist2))
            hb_dist[k] = -hb_dist[k];
    }

    /* post processing to account for cases such as (I34 T--A J259) 1zla -- nov-02-2006 */
    for (m = 1; m <= num_hbonds; m++) {
        matched_idx[m] = 0;
        k = idx2[m][1] + idx2[m][2];
        if (k != 18)
            continue;
        for (n = 1; n <= num_hbonds; n++) {
            if (hb_dist[n] < 0.0)  /* covering n == m, 9 + 9 == 18, 0 + 1 etc */
                continue;
            if (hb_dist[n] > fabs(hb_dist[m]) ||
                !dval_in_range(hb_dist[n], misc_pars->hb_lower, misc_pars->hb_dist2))
                continue;
            if (!strcmp(hb_atom1[m], hb_atom1[n]) || !strcmp(hb_atom2[m], hb_atom2[n])) {
                hb_dist[n] = -hb_dist[n];
                matched_idx[n] = 1;
            }
        }
    }

    for (k = 1; k <= num_hbonds; k++)
        if (matched_idx[k] && !dval_in_range(fabs(hb_dist[k]), misc_pars->hb_lower,
                                             misc_pars->hb_dist2))
            hb_dist[k] = fabs(hb_dist[k]);
#endif
/* ============================================================================ */

    free_lmatrix(idx2, 1, num_hbonds, 1, 2);
    free_lvector(matched_idx, 1, num_hbonds);

    /* === End of further processing by adding more possible H-bonds === */
}

long validate_hbonds(long num_hbonds, double *hb_dist, long *lkg_type, char *hb_type,
                     char basei, char basej, char **hb_atom1, char **hb_atom2)
{
    long k, m = 0;

    for (k = 1; k <= num_hbonds; k++) {
        hb_type[k] = ' ';  /* initialized to ' ' */
        if (hb_dist[k] > 0.0)  /* not taken as H-bond */
            continue;
        hb_type[k] = donor_acceptor(basei, basej, hb_atom1[k], hb_atom2[k]);
        hb_dist[k] = fabs(hb_dist[k]);  /* using hb_type as an indicator */
        if (hb_type[k] == '-' && dval_in_range(hb_dist[k], 2.5, 3.5))
            m++;
    }

    if (m) {  /* with good D-A type H-bond, further checking & elimination */
        for (k = 1; k <= num_hbonds; k++) {
            if (hb_type[k] == ' ')
                continue;
            if ((hb_dist[k] > 3.6) ||  /* must be within reasonable limit 3.6 */
                (hb_type[k] == '*' && lkg_type[k] != 18 &&  /* D-D or A-A, not in best link */
                 !dval_in_range(hb_dist[k], 2.6, 3.2)))  /* e.g. A1498[U-**+-A]A1499 of 1vs7 */
                hb_type[k] = ' ';
        }
    }

    /* now count and output 'good' H-bonding information */
    m = 0;
    for (k = 1; k <= num_hbonds; k++) {
        if (hb_type[k] == ' ')
            continue;
        m++;
    }

    return m;
}

/* get H-bond length information between residue i and j */
void get_hbond_ij(long i, long j, char basei, char basej, miscPars * misc_pars,
                  long **seidx, long *idx, char **AtomName, double **xyz, char *hb_info)
{
    char *hb_type, **hb_atom1, **hb_atom2, aname1[5], aname2[5], stmp[20];
    double *hb_dist;
    long k, m, n, num_hbonds = 0, *lkg_type;

    hb_atom1 = cmatrix(1, BUF512, 0, 4);
    hb_atom2 = cmatrix(1, BUF512, 0, 4);
    hb_dist = dvector(1, BUF512);

    for (m = seidx[i][1]; m <= seidx[i][2]; m++) {
        for (n = seidx[j][1]; n <= seidx[j][2]; n++) {
            if (good_hbatoms(misc_pars, AtomName[m], AtomName[n], idx[m], idx[n]) &&
                within_limits(xyz[n], xyz[m], misc_pars->hb_lower, misc_pars->hb_dist1)) {
                if (++num_hbonds > BUF512)
                    fatal("Too many possible H-bonds between two bases\n");
                strcpy(hb_atom1[num_hbonds], AtomName[m]);
                strcpy(hb_atom2[num_hbonds], AtomName[n]);
                hb_dist[num_hbonds] = p1p2_dist(xyz[n], xyz[m]);
            }
        }
    }

    /* Yurong reported a case in rr0027 (1i97): >A:..39_:[..G]G-**+-A[..A]:.530_:A<
       where the bp criteria defines it as a pair, because G:N7-A:N7 distance is 3.89A
       which is within 4.0A cut-off. However, in this function when H-bonding info is
       checked, there is no good_hbatoms() pairs. Thus num_hbonds is ZERO and
       hb_info[] would be undefined, leaving it to whatever content the system assign
       it with -- mostly garbage. Here we explicitly initialize it to an empty string to
       avoid this problem. */
    if (!num_hbonds)  /* no H-bond found */
        sprintf(hb_info, "[%ld]", num_hbonds);

    else {
        lkg_type = lvector(1, num_hbonds);
        hb_atompair(num_hbonds, hb_atom1, hb_atom2, hb_dist, lkg_type, misc_pars);

        /*  to eliminate erroneous H-bonds, based purely on geometrical criteria */
        hb_type = cvector(1, num_hbonds);
        m = validate_hbonds(num_hbonds, hb_dist, lkg_type, hb_type, basei, basej,
                            hb_atom1, hb_atom2);

        sprintf(hb_info, "[%ld]", m);
        for (k = 1; k <= num_hbonds; k++) {
            if (hb_type[k] == ' ')
                continue;
            strcpy(aname1, hb_atom1[k]);
            cvt_pdbv3_name(aname1);
            strcpy(aname2, hb_atom2[k]);
            cvt_pdbv3_name(aname2);
            sprintf(stmp, " %s%c%s %4.2f", aname1, hb_type[k], aname2, hb_dist[k]);
            strcat(hb_info, stmp);
        }

        free_cvector(hb_type, 1, num_hbonds);
        free_lvector(lkg_type, 1, num_hbonds);
    }

    free_cmatrix(hb_atom1, 1, BUF512, 0, 4);
    free_cmatrix(hb_atom2, 1, BUF512, 0, 4);
    free_dvector(hb_dist, 1, BUF512);
}

/* based on canonical form of ACGITU, decide if a HB-pair is donor-acceptor */
char donor_acceptor(char basei, char basej, char *hb_atom1, char *hb_atom2)
{
    char da[3], hbatom_type = '*', *pchar;
    char ia = '\0', ja = '\0';
    /* _A: acceptor; _D: donor; _X: both acceptor & donor; _?: not sure */
    static char *cmn_base = CB_LIST;
    static char *da_type[7] = { "AD", "AX", "XD", "XX", "DA", "DX", "XA" };
    static char *bb_da[6] = { " O1P_A", " O2P_A", " O5'_A", " O4'_A", " O3'_A",
        " O2'_X"
    };
    static char *base_da[6][6] = {
        {" N9 _?", " N7 _A", " N6 _D", " N1 _A", " N3 _A"},  /* A */
        {" N1 _?", " O2 _A", " N3 _A", " N4 _D"},  /* C */
        {" N9 _?", " N7 _A", " O6 _A", " N1 _D", " N2 _D", " N3 _A"},  /* G */
        {" N9 _?", " N7 _A", " O6 _A", " N1 _D", " N3 _A"},  /* I */
        {" N1 _?", " O2 _A", " N3 _D", " O4 _A"},  /* T */
        {" N1 _?", " O2 _A", " N3 _D", " O4 _A"}  /* U */
    };
    long i, inum = -1, jnum = -1, num = 6;

    if ((pchar = strchr(cmn_base, toupper((int) basei))) != NULL)
        inum = pchar - cmn_base;
    if ((pchar = strchr(cmn_base, toupper((int) basej))) != NULL)
        jnum = pchar - cmn_base;
    if (inum >= 0 && jnum >= 0) {  /* both are one of the six bases */
        for (i = 0; i < num; i++) {  /* check if backbone atom */
            if (!strncmp(bb_da[i], hb_atom1, 4))
                ia = bb_da[i][5];
            if (!strncmp(bb_da[i], hb_atom2, 4))
                ja = bb_da[i][5];
        }
        if (!ia)
            for (i = 0; i < num; i++)
                if (base_da[inum][i] && !strncmp(base_da[inum][i], hb_atom1, 4)) {
                    ia = *(base_da[inum][i] + 5);
                    break;
                }
        if (!ja)
            for (i = 0; i < num; i++)
                if (base_da[jnum][i] && !strncmp(base_da[jnum][i], hb_atom2, 4)) {
                    ja = *(base_da[jnum][i] + 5);
                    break;
                }
        if (ia && ja) {
            sprintf(da, "%c%c", ia, ja);
            if (num_strmatch(da, da_type, 0, 6))
                hbatom_type = '-';
        }
    }

    return hbatom_type;
}

long asym_idx(char *asym, char atoms_list[NELE][3], long dft_lval)
{
    long i;

    for (i = 0; i < NELE; i++)
        if (!strcmp(asym, atoms_list[i]))
            return i;

    return dft_lval;
}

void atom_info(long idx, char atoms_list[NELE][3], double *covalence_radii, double *vdw_radii)
{
    static char *ALIST[NELE] =
        { UNKATM, " C", " O", " H", " N", " S", " P", " F", "CL", "BR", " I",
        "SI"
    };
    static double CRADII[NELE] = { 1.666, 0.762, 0.646, 0.352, 0.689, 1.105, 1.000, 0.619, 1.022,
        1.183, 1.378, 1.105  /* XX original 1.200 */
    };
    static double VRADII[NELE] =
        { 2.00, 1.70, 1.52, 1.20, 1.55, 1.80, 1.80, 1.47, 1.75, 1.85, 1.98,
        2.10
    };
    long i;

    if (idx == 1)
        for (i = 0; i < NELE; i++)
            strcpy(atoms_list[i], ALIST[i]);
    else if (idx == 2)
        for (i = 0; i < NELE; i++)
            covalence_radii[i] = CRADII[i];
    else if (idx == 3)
        for (i = 0; i < NELE; i++)
            vdw_radii[i] = VRADII[i];
    else
        fatal("wrong options for <atom_info>: should be 1, 2 or 3\n");
}

/* get 2-letter atomic symbol from atom name in PDB file */
void aname2asym(const char *aname0, char *my_asym, long num_sa, char **atomlist)
{
    char aname[BUF512];
    long i, unknown, NAME_LEN = 4;

    strcpy(aname, aname0);

    /* non-alphabets changed to '.': e.g., " O1P": ".O.P"; " N1 ": ".N.." */
    for (i = 0; i < NAME_LEN; i++)
        if (!isalpha((int) aname[i]))
            aname[i] = '.';

    for (i = 1; i <= num_sa; i++)
        if (str_pmatch(atomlist[i], aname))
            break;

    if (i > num_sa) {
        unknown = is_equal_string(aname, ".UNK");

        if ((aname[0] != '.') && (aname[1] != '.') && (aname[2] == '.') && (aname[3] == '.')) {  /* as in 'NA..' */
            my_asym[0] = aname[0];
            my_asym[1] = aname[1];
            my_asym[2] = '\0';
        } else if ((aname[0] == '.') && (aname[1] != '.') && !unknown) {
            my_asym[0] = ' ';
            my_asym[1] = aname[1];
            my_asym[2] = '\0';
        } else if (aname[0] == 'H')
            strcpy(my_asym, " H");
        else
            strcpy(my_asym, UNKATM);

        if (!unknown) {
            fprintf(stderr, "no matching entry for atom name [%s] (%s) in '%s'\n",
                    aname0, aname, ATOM_FILE);
            fprintf(stderr, "\tnow it is set as '%s'\n", my_asym);
            fprintf(stderr, "\tcheck and update file $X3DNA/config/atomlist.dat\n");
        }
    } else
        strcpy(my_asym, atomlist[i] + NAME_LEN);
}

/* get atom index for calculating bond linkage */
void atom_idx(long num, char **AtomName, char **Miscs, long *idx)
{
    char pdb_asym[3], my_asym[3], atoms_list[NELE][3];
    long i, k, bad_pdb_asym = FALSE;

    atom_info(1, atoms_list, NULL, NULL);

    for (i = 1; i <= num; i++) {
        strcpy(pdb_asym, UNKATM);  /* atomic symbol from PDB data file */
        k = FALSE;  /* default to no atomic symbol */

        if (Miscs && strlen(Miscs[i]) >= 27 && !str_pmatch(Miscs[i] + 25, "  ")) {
            strncpy(pdb_asym, Miscs[i] + 25, 2);  /* with this info. */
            pdb_asym[2] = '\0';

            if (is_equal_string(pdb_asym, " D"))
                strcpy(pdb_asym, " H");

            if (num_strmatch(pdb_asym, Gvars.ATOM_NAMES, 0, Gvars.NUM_ELE))
                k = TRUE;
            else
                bad_pdb_asym = TRUE;  /* for overall checking */
        }

        if (k && !bad_pdb_asym)
            strcpy(my_asym, pdb_asym);
        else  /* atomic symbol deduced from atom name */
            aname2asym(AtomName[i], my_asym, Gvars.NUM_SATOM, Gvars.ATOMLIST);

        idx[i] = asym_idx(my_asym, atoms_list, 0);

        /* fprintf(stderr, "%ld\t%s\t%s\t%ld\n", i, AtomName[i], my_asym, idx[i]); */
    }

    if (bad_pdb_asym)
        fprintf(stderr, "\tPDB with illegal atomic symbol in columns #77-78\n");
}

/* get linkage information based on inter-atomic distance */
void get_bonds(long num, char **AtomName, double **xyz, long num_residue, long *RY,
               long **seidx, long **connect)
{
    long i, nbond, nbond_estimated, *idx, **linkage;

    idx = lvector(1, num);
    atom_idx(num, AtomName, NULL, idx);

    nbond_estimated = lround(NBOND_FNUM * NUM_RESIDUE_ATOMS);
    linkage = lmatrix(1, nbond_estimated, 1, 2);
    for (i = 1; i <= num_residue; i++)  /* linkage within each residue */
        if (RY[i] >= 0) {  /* a base residue */
            nbond = 0;
            atom_linkage(seidx[i][1], seidx[i][2], idx, xyz, NULL, NULL,
                         nbond_estimated, &nbond, linkage);
            lkg2connect(AtomName, seidx[i][1], seidx[i][2], nbond, linkage, connect);
        }

    free_lvector(idx, 1, num);
    free_lmatrix(linkage, 1, nbond_estimated, 1, 2);
}

/* get atom linkage using covalent radii criterion. alt. position & chain info. */
void atom_linkage(long ib, long ie, long *idx, double **xyz, char **Miscs, char *ChainID,
                  long nbond_estimated, long *nbond, long **linkage)
{
    char a1, a2, c1, c2;
    double dst, covalence_radii[NELE];
    long j, k;

    /* Bond criteria: 1.15 * (rA + rB)
       RasMol: 0.56 + (rA + rB)   MacroModel: 1.25 * (rA + rB) */

    atom_info(2, NULL, covalence_radii, NULL);
    for (j = ib; j <= ie - 1; j++) {
        a1 = (Miscs == NULL) ? ' ' : Miscs[j][1];  /* atom alternative position */
        c1 = (ChainID == NULL) ? '_' : ChainID[j];
        for (k = j + 1; k <= ie; k++) {
            a2 = (Miscs == NULL) ? ' ' : Miscs[k][1];  /* atom alternative position */
            c2 = (ChainID == NULL) ? '_' : ChainID[k];
            if (a1 != ' ' && a2 != ' ' && a1 != a2)  /* same model: alt. position */
                continue;
            if (c1 != c2)  /* same chain ID */
                continue;
            dst = BOND_FACTOR * (covalence_radii[idx[j]] + covalence_radii[idx[k]]);
            if (within_limits(xyz[k], xyz[j], 0, dst)) {
                if (++*nbond > nbond_estimated)
                    fatal("too many linkages\n");
                else {
                    linkage[*nbond][1] = j;
                    linkage[*nbond][2] = k;
                }
            }
        }
    }
}

/* from bond linkage to PDB connection table: SLOW for large structure rr0033 */
void lkg2connect(char **AtomName, long ib, long ie, long nbond, long **linkage, long **connect)
{
    long i, ilink, j;
    long idx[7];

    for (i = ib; i <= ie; i++) {
        ilink = 0;
        for (j = 1; j <= nbond; j++) {
            if (i == linkage[j][1]) {
                if (++ilink > 6) {
                    fprintf(stderr, "atom <%ld: [%s]> has over 6 bonds\n", i, AtomName[i]);
                    break;
                } else
                    idx[ilink] = linkage[j][2];
            }
            if (i == linkage[j][2]) {
                if (++ilink > 6) {
                    fprintf(stderr, "atom <%ld: [%s]> has over 6 bonds\n", i, AtomName[i]);
                    break;
                } else
                    idx[ilink] = linkage[j][1];
            }
        }
        if (ilink > 6)
            ilink = 6;  /* maximum six bonds */
        if (ilink > 0) {
            for (j = 1; j <= ilink; j++)
                connect[i][j] = idx[j];
            connect[i][7] = ilink;
            lsort(ilink, connect[i], idx);  /* ignore idx */
        }
    }
}

/* initialize matrix htm_water as follows:
 * row#1: c0 for # of atoms; c1-num for sequential residue # in the PDB file
 * row#2: c0 for # of residues; c1-num for connected base sequential residue #
 * row#3: c0 for hydration indication; c1-num for atom index
 * row#4: c0 for # of waters; c1--num_H2O for corresponding residue # */
void init_htm_water(long waters, long num, long num_residue, long *idx, long **htm_water)
{
    long i;

    htm_water[1][0] = num;  /* total # of atoms */
    htm_water[2][0] = num_residue;  /* total # of residue */
    htm_water[3][0] = waters;  /* if to check for hydration */

    for (i = 1; i <= num; i++)  /* make a copy of atom index */
        htm_water[3][i] = idx[i];
}

/* make each HETATM record, including water, to its connected base residue */
void identify_htw(long num_residue, long **seidx, long *RY, char **AtomName,
                  char **ResName, char *ChainID, long *ResSeq, char **Miscs,
                  double **xyz, long **htm_water)
{
    static char *WATER[] = { WATER_LIST };
    char a1, a2, c1, c2;
    double dst, covalence_radii[NELE];
    long num_wat = sizeof WATER / sizeof WATER[0] - 1;
    long i, ib, ie, id, j, k = 0, m, num_H2O = 0;

    /* identify ligands and water molecules */
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] >= 0)  /* skip base residues */
            continue;
        ib = seidx[i][1];
        ie = seidx[i][2];
        id = ie - ib + 1;
        for (j = ib; j <= ie; j++) {
            if (Miscs[j][0] == 'H') {  /* HETATM record */
                htm_water[1][j] = i;
                htm_water[2][j] = -8;  /* general HETATM record indicator */
            }
        }

        if (num_strmatch(ResName[ib], WATER, 0, num_wat)  /* matching H2O ResName */
            ||(id == 1 &&  /* ONLY ONE atom in this residue */
               ((htm_water[3][ib] == 2)  /* added nov-03-2006: 2 is Oxygen */
                ||!strcmp(AtomName[ib], " O  ") || !strcmp(AtomName[ib], " OW ")))) {
            for (j = ib; j <= ie; j++) {  /* possibly with Hs, e.g. NMR */
                htm_water[1][j] = i;
                htm_water[2][j] = -1;  /* for water molecule */
            }
            num_H2O++;
            htm_water[4][num_H2O] = i;
        }
    }
    htm_water[4][0] = num_H2O;

    /* This section checks for covalent connection of a non-base residue,
     * e.g., drug molecule, to a base residue. They have to be on the same
     * alternative position and chain ID. Pure H2Os do not count here. */
    atom_info(2, NULL, covalence_radii, NULL);  /* cf. <atom_linkage> */
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] >= 0)  /* skip base residues */
            continue;
        id = 0;  /* check if residue "i" is connected to a base residue */
        for (j = seidx[i][1]; j <= seidx[i][2]; j++) {
            /* assuming NO more than 1 connection to base; only HETATM/H2O */
            if (htm_water[2][j] >= 0)  /* -8 or -1 */
                break;
            a1 = Miscs[j][1];
            c1 = ChainID[j];
            for (k = 1; k <= num_residue; k++) {
                if (RY[k] < 0)  /* check only for connections to base residue */
                    continue;
                for (m = seidx[k][1]; m <= seidx[k][2]; m++) {
                    a2 = Miscs[m][1];
                    c2 = ChainID[m];
                    if (c1 != c2)  /* same chain ID */
                        break;  /* restricted to residue m */
                    if (a1 != ' ' && a2 != ' ' && a1 != a2)
                        continue;  /* same model: alt. position */
                    dst = BOND_FACTOR * (covalence_radii[htm_water[3][j]] +
                                         covalence_radii[htm_water[3][m]]);
                    if (within_limits(xyz[m], xyz[j], 0, dst)) {
                        id = 1;
                        goto ID_CNCT;
                    }
                }
            }
        }

      ID_CNCT:
        if (id) {  /* residue i connected to base residue k */
            for (j = seidx[i][1]; j <= seidx[i][2]; j++)
                htm_water[2][j] = k;
            ib = seidx[i][1];  /* HETATM residue */
            ie = seidx[k][1];  /* contacted base residue */

            /* no longer check for match of residue #: ResSeq[ib] != ResSeq[ie]: 4tna */
            if (Miscs[ib][2] != Miscs[ie][2])  /* e.g., 1fjb: 17(697,23) and 17(507,17) */
                fprintf(stderr, "%ld(%ld,%ld) and %ld(%ld,%ld) have different"
                        " insertion code!\n", ResSeq[ib], ib, i, ResSeq[ie], ie, k);
        }
    }
}

/* attach covalently connected residues as in +A etc */
long attached_residues(long inum_base, long *ivec, long *ivec2, long **seidx,
                       double **xyz, long **htm_water, miscPars * misc_pars)
{
    long i, iw, ir, j, k, m, n, tnum_res = inum_base;
    long num_residue = htm_water[2][0], num_H2O = htm_water[4][0];

    /* check if there is a covalent bond linkage between hetero
     * group (e.g., a drug molecule) with the base residue. */
    for (i = 1; i <= inum_base; i++) {
        ivec2[i] = ivec[i];

        if (Gvars.ATTACH_RESIDUE == FALSE)  /* Per Pascal's request: no metal, or HETATM */
            continue;

        for (j = 1; j <= num_residue; j++) {
            k = seidx[j][1];
            if (htm_water[2][k] == ivec[i])
                ivec2[++tnum_res] = htm_water[1][k];
        }
    }

    if (!htm_water[3][0])
        return tnum_res;

    inum_base = tnum_res;
    n = inum_base + 1;  /* starting index for hydration residues */
    for (ir = 1; ir <= num_H2O; ir++) {  /* time-consuming part: best way! */
        iw = htm_water[4][ir];  /* hydration index */
        for (m = seidx[iw][1]; m <= seidx[iw][2]; m++) {
            if (htm_water[3][m] != 2)  /* not O */
                continue;
            for (i = 1; i <= inum_base; i++) {
                k = ivec2[i];
                for (j = seidx[k][1]; j <= seidx[k][2]; j++) {
                    if (!lval_in_set(htm_water[3][j], 1, misc_pars->water_idx[0],
                                     misc_pars->water_idx))
                        continue;
                    if (within_limits(xyz[m], xyz[j], misc_pars->hb_lower, misc_pars->water_dist)
                        && !lval_in_set(iw, n, tnum_res, ivec2)) {
                        ivec2[++tnum_res] = iw;
                    }
                }
            }
        }
    }

    if (tnum_res > inum_base + 1) {  /* sort the hydration residues into order */
        long *idx;

        k = tnum_res - inum_base;
        idx = lvector(1, k);
        lsort(k, ivec2 + inum_base, idx);

        free_lvector(idx, 1, k);
    }

    return tnum_res;
}

/* print out base-pairing information */
void print_pairinfo(long i, long j, char basei, char basej, double *rtn_val, double *chi,
                    miscPars * misc_pars, long **seidx, long *idx, char **AtomName,
                    double **xyz, char *bseq, long detailed, FILE * fp)
{
    char antip, bptype[4], hb_info[BUF512];
    long k;

    antip = (rtn_val[35] < 0.0) ? '-' : '+';
    sprintf(bptype, "%c%c%c", bseq[i], antip, bseq[j]);
    fprintf(fp, "             %s  ", bptype);
    for (k = 27; k <= 32; k++)  /* bp parameters */
        fprintf(fp, "%8.2f", rtn_val[k]);
    fprintf(fp, "\n");
    get_hbond_ij(i, j, basei, basej, misc_pars, seidx, idx, AtomName, xyz, hb_info);
    fprintf(fp, "             %s %s\n", bptype, hb_info);

    if (!detailed)
        return;

    fprintf(fp, "                    %s %s %s\n", rtn_val[35] < 0.0 ? "anti-parallel" : "parallel",  /* z-axis */
            rtn_val[33] < 0.0 ? "trans" : "cis",  /* x-axis */
            rtn_val[36] < 0.0 ? "trans" : "cis");  /* C1-N vector */
    fprintf(fp, "            ");
    for (k = 33; k <= 35; k++)
        fprintf(fp, "%c", rtn_val[k] > 0.0 ? '+' : '-');
    for (k = 33; k <= 36; k++)
        fprintf(fp, "%9.1f", dot2ang(rtn_val[k]));
    fprintf(fp, "%9.1f%9.1f\n", chi[i], chi[j]);

    fprintf(fp, "     ");
    for (k = 1; k <= 8; k++)
        fprintf(fp, "%8.2f", rtn_val[k]);
    fprintf(fp, "\n");
}

/* once a pair is found, calculate more parameters */
static void calculate_more_bppars(long i, long j, double dir_x, double dir_y,
                                  double dir_z, double **orien, double **org, char *bseq,
                                  double **NC1xyz, double *rtn_val, long *bpid)
{
    char bpi[3];
    long k, l, koffset;
    double zave[4], dNN_vec[4], pars[7], **r1, **r2, **mst;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);

    for (k = 1; k <= 9; k++) {
        rtn_val[k + 8] = orien[i][k];  /* base I xyz axes */
        rtn_val[k + 17] = orien[j][k];  /* base II xyz axes */
    }
    for (k = 1; k <= 3; k++) {
        koffset = (k - 1) * 3;
        for (l = 1; l <= 3; l++) {
            r1[l][k] = orien[i][koffset + l];
            r2[l][k] = (k == 1 || dir_z > 0) ?  /* keep x, reverse y & z if anti-parallel */
                orien[j][koffset + l] : -orien[j][koffset + l];
        }
    }
    bpstep_par(r2, org[j], r1, org[i], pars, mst, &rtn_val[5]);
    for (k = 1; k <= 6; k++)  /* bp parameters in columns 27-32 */
        rtn_val[26 + k] = pars[k];

    sprintf(bpi, "%c%c", toupper((int) bseq[i]), toupper((int) bseq[j]));
    *bpid = -1;  /* assumed to be non-WC */
    if (dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0) {
        check_wc_wobble_pair(bpid, bpi, pars[1], pars[2], pars[6]);
        if (*bpid == 2)
            rtn_val[5] -= 2.0;  /* WAS: -1.5 (bonus for WC pair) */
    }

    rtn_val[33] = dir_x;
    rtn_val[34] = dir_y;
    rtn_val[35] = dir_z;

    if (NC1xyz[i][7] > 0 && NC1xyz[j][7] > 0) {
        ddxyz(NC1xyz[i], NC1xyz[i] + 3, zave);
        ddxyz(NC1xyz[j], NC1xyz[j] + 3, dNN_vec);
        vec_norm(zave);
        vec_norm(dNN_vec);
        rtn_val[36] = dot(zave, dNN_vec);
    } else
        rtn_val[36] = EMPTY_NUMBER;

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
}

/* for judging the quality of a pair. mostly used to find the "best" pair */
static double adjust_pairQuality(long i, long j, char basei, char basej, long **seidx,
                                 long *idx, char **AtomName, double **xyz, miscPars * misc_pars)
{
    double dval;
    long k, num_hb, num_good_hb = 0, **num_list;

    num_list = lmatrix(1, BUF512, 0, 3);
    hb_numlist(i, j, basei, basej, seidx, idx, AtomName, xyz, misc_pars, &num_hb, num_list);
    for (k = 1; k <= num_hb; k++) {
        if (num_list[k][0])  /* not canonical D-A type */
            continue;
        dval = num_list[k][3] / MFACTOR;
        if (dval_in_range(dval, 2.5, 3.5))  /* assumed good H-bonding distance */
            num_good_hb++;
    }
    free_lmatrix(num_list, 1, BUF512, 0, 3);

    if (num_good_hb >= 2)  /* make it works for i54:55 vs j238:239 of 1zla */
        return -3.0;
    else
        return -num_good_hb;
}

/* Checking if two bases form a pair according to several criteria
 * rtn_val[RTNNUM]:
 *       d, dv, angle, dNN, dsum, bp-org,  x1,  y1,   z1,   x2,   y2,   z2
 * col#  1   2    3     4     5    6-8    9-11 12-14 15-17 18-20 21-23 24-26
 * bpid: 0: not-paired; +1: WC geometry; +2: WC pair; -1: other cases
 *       bp-pars,   bp_relative_orientation, C1_N_relative_orientation
 * col#   27-32        33-35 (dx, dy, dz)              36
 *
 * There are also possibilities when the pair is mostly maintained by
 * backbone + backbone; backbone + base H-bonds. e.g. A512 + C637
 * in rr0033 (chain 0). No direct base-to-base H-bonds ************** */
void check_pair(long i, long j, char *bseq, long **seidx, double **xyz,
                double **NC1xyz, double **orien, double **org, long *idx,
                char **AtomName, miscPars * misc_pars, double *rtn_val,
                long *bpid, long **ring_atom, long network)
{
    double dir_x, dir_y, dir_z;
    double dorg[4], oave[4], zave[4], dNN_vec[4];
    long cdns, m, n, num_base_hb = 0, num_o2_hb = 0;

    *bpid = 0;  /* default as not-paired */
    if (i == j)
        return;  /* same residue */

    get_bp_zoave(i, j, orien, org, oave, zave);

    ddxyz(org[i], org[j], dorg);
    ddxyz(NC1xyz[i], NC1xyz[j], dNN_vec);

    rtn_val[1] = veclen(dorg);  /* distance between origins */
    dir_x = dot(&orien[i][0], &orien[j][0]);  /* relative x direction */
    dir_y = dot(&orien[i][3], &orien[j][3]);  /* relative y direction */
    dir_z = dot(&orien[i][6], &orien[j][6]);  /* relative z direction */

    rtn_val[2] = fabs(dot(dorg, zave));  /* dv: projection onto mean normal */
    rtn_val[3] = z1_z2_angle_in_0_to_90(&orien[i][6], &orien[j][6]);  /* angle between base normals */
    rtn_val[4] = veclen(dNN_vec);  /* RN9-YN1 distance */
    rtn_val[5] = rtn_val[1] + 2.0 * rtn_val[2] + rtn_val[3] / 20.0;
    /* WAS 25. changed based on 1id3 i25:26/j267:268; i107:8/j185:6 */

    if (network) {  /* check if two bases in pairing network */
#if 0
        fprintf(stderr, "i: %ld\tj: %ld\tangle: %g\n", i, j, rtn_val[3]);
#endif
        if (dval_in_range(rtn_val[3], misc_pars->min_plane_angle, misc_pars->max_plane_angle)
            && dval_in_range(rtn_val[4], misc_pars->min_dNN, misc_pars->max_dNN)
/* -- more strict than for best-pair below; otherwise, np_recipes/R5_pentaplets
      will have different result at least for 'multiplets.pdb'
            && get_oarea(i, j, ring_atom, oave, zave, xyz, 0) < OVERLAP)
*/
            && get_oarea(i, j, ring_atom, oave, zave, xyz, 0) < OVERLAP &&
            (get_oarea(i, j, ring_atom, org[i], &orien[i][6], xyz, 0) < OVERLAP) &&
            (get_oarea(i, j, ring_atom, org[j], &orien[j][6], xyz, 0) < OVERLAP))
            *bpid = -1;  /* WAS 1, Oct. 30, 2006 */
        rtn_val[5] += adjust_pairQuality(i, j, bseq[i], bseq[j], seidx, idx, AtomName,
                                         xyz, misc_pars);
        return;
    }

    cdns = (dval_in_range(rtn_val[1], misc_pars->min_dorg, misc_pars->max_dorg) &&
            dval_in_range(rtn_val[2], misc_pars->min_dv, misc_pars->max_dv) &&
            dval_in_range(rtn_val[3], misc_pars->min_plane_angle, misc_pars->max_plane_angle)
            && dval_in_range(rtn_val[4], misc_pars->min_dNN, misc_pars->max_dNN));
/* no-overlap with reference to middle, lower (i) or upper (j) frame
    if (cdns &&
        (get_oarea(i, j, ring_atom, oave, zave, xyz, 0) < OVERLAP) &&
        (get_oarea(i, j, ring_atom, org[i], &orien[i][6], xyz, 0) < OVERLAP) &&
        (get_oarea(i, j, ring_atom, org[j], &orien[j][6], xyz, 0) < OVERLAP)) {
*/
    if (cdns && get_oarea(i, j, ring_atom, oave, zave, xyz, 0) < OVERLAP) {
        for (m = seidx[i][1]; m <= seidx[i][2]; m++)
            for (n = seidx[j][1]; n <= seidx[j][2]; n++) {
                if (!within_limits(xyz[m], xyz[n], misc_pars->hb_lower, misc_pars->hb_dist1))
                    continue;

                if (is_baseatom(AtomName[m]) && is_baseatom(AtomName[n]) &&
                    good_hbatoms(misc_pars, AtomName[m], AtomName[n], idx[m], idx[n]))
                    num_base_hb++;

                if (is_equal_string(AtomName[m], " O2'") || is_equal_string(AtomName[n], " O2'"))
                    num_o2_hb++;
            }

        if ((misc_pars->min_base_hb && (num_base_hb >= misc_pars->min_base_hb)) ||
            (!misc_pars->min_base_hb && (num_o2_hb || num_base_hb))) {
            calculate_more_bppars(i, j, dir_x, dir_y, dir_z, orien, org, bseq, NC1xyz,
                                  rtn_val, bpid);
            rtn_val[5] +=
                adjust_pairQuality(i, j, bseq[i], bseq[j], seidx, idx, AtomName, xyz, misc_pars);
        }
    }
}

/* get the xyz coordinates of O3' and P atoms */
void o3_p_xyz(long ib, long ie, char *aname, char **AtomName, double **xyz,
              double *o3_or_p, long idx)
{
    long i;

    i = find_1st_atom(aname, AtomName, ib, ie, "");
    if (i) {  /* with O3'/P atom */
        cpxyz(xyz[i], o3_or_p + idx - 4);
        o3_or_p[idx] = 1.0;
    } else  /* without O3'/P atom */
        o3_or_p[idx] = -1.0;
}

/* check if atomname is a nucleotide base atom: ASSUMING NA residue */
long is_baseatom(char *atomname)
{
    if (is_equal_string(atomname, " C5M"))  /* C5M of T */
        return TRUE;

    if (atomname[0] == ' ' && strchr("HP", atomname[1]) == NULL  /* like " N1 " */
        && isdigit((int) atomname[2]) && atomname[3] == ' ')
        return TRUE;

    return FALSE;
}

static long glyco_N(long isR, long ib, long ie, char b, char **AtomName, char **ResName,
                    long C1prime, double **xyz)
{
    char *a;
    double d, dm = 9999.0, *c1xyz;
    long k, km, num = 0;

    if (isR) {  /* purine */
        k = find_1st_atom(" N9 ", AtomName, ib, ie, "");
        if (k)
            return k;
    } else {  /* pyrimidine */
        if (strchr("Pp", b)) {  /* pseudo-uridine */
            k = find_1st_atom(" C5 ", AtomName, ib, ie, "");
            if (k)
                return k;
        } else {  /* regular pyrimidine */
            k = find_1st_atom(" N1 ", AtomName, ib, ie, "");
            if (k)
                return k;
        }
    }

    assert(!k);
    fprintf(stderr, "Cannot identify RN9/YN1 in [%s] -- ", ResName[ib]);

    if (C1prime) {  /* find the shorest distance */
        c1xyz = xyz[C1prime];
        km = FALSE;
        for (k = ib; k <= ie; k++) {
            a = AtomName[k];
            if (!is_baseatom(a))
                continue;
            d = p1p2_dist(c1xyz, xyz[k]);
            if (d < dm) {
                dm = d;
                km = k;
            }
        }
    }

    if (dm <= BOND_UPPER_LIMIT) {
        fprintf(stderr, "use atom [%s] instead\n", AtomName[km]);
        return km;
    }

    /* last try! */
    km = FALSE;
    for (k = ib; k <= ie; k++) {
        a = AtomName[k];
        if ((isR && strchr(a, '9')) || (!isR && strchr(a, '1'))) {
            num++;
            km = k;
        }
    }

    if (num == 1) {
        fprintf(stderr, "[i] using atom [%s] instead\n", AtomName[km]);
        return km;
    }

    fatal("stop!\n");

    return 0;
}

/* get base information for locating possible pairs later */
void base_info(long num_residue, char *bseq, long **seidx, long *RY,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, double **orien, double **org,
               double **NC1xyz, double **o3_p)
{
    char BDIR[BUF512];
    long i, ib, ie, C1prime, N;

    get_BDIR(BDIR, "Atomic_A.pdb");

    /* get the reference frame for each base */
    base_frame(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
               ResSeq, Miscs, xyz, BDIR, orien, org);

    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];
        if (RY[i] >= 0) {  /* a base */
            C1prime = find_1st_atom(" C1'", AtomName, ib, ie, "");
            if (C1prime) {
                NC1xyz[i][7] = 1.0;  /* C1' atom exists */
                cpxyz(xyz[C1prime], NC1xyz[i] + 3);  /* 4-6 */
            }
            N = glyco_N(RY[i], ib, ie, bseq[i], AtomName, ResName, C1prime, xyz);
            cpxyz(xyz[N], NC1xyz[i]);  /* RN9/YN1 surely exist */
            o3_p_xyz(ib, ie, " O3'", AtomName, xyz, o3_p[i], 4);
            o3_p_xyz(ib, ie, " P  ", AtomName, xyz, o3_p[i], 8);
        }
    }
}

void help3dna_usage(char *program_name)
{
    help3dna(program_name);
    contact_msg(0);
}

/* Easy 3DNA help message read from a text file */
void help3dna(char *program_name)
{
    char BDIR[BUF512], str[BUF512], pname[BUF512], *prefix = "        ";
    long nlen, nchar = 75, found_help = 0;
    FILE *fp;

    get_BDIR(BDIR, HELP3DNA);
    strcat(BDIR, HELP3DNA);
    fp = open_file(BDIR, "r");

    print_sep(stderr, '=', nchar);
    strcpy(pname, program_name);  /* program_name is a constant */
    nlen = upperstr(pname);  /* to upper case */
    while (fgets(str, sizeof str, fp) != NULL) {
        if (str[0] != '<')
            continue;
        upperstr(str);
        if (strncmp(str + 1, pname, nlen))
            continue;  /* not match */
        found_help = 1;
        while (fgets(str, sizeof str, fp) != NULL) {
            if (str[0] == '<') {
                upperstr(str);
                if (str[1] != '/' || (strncmp(str + 2, pname, nlen)))
                    fatal("error in help format: <tag> ... </tag>\n");
                found_help = 2;
                goto FINISHED;
            } else {
                if (str[0] != '#')  /* not-comment */
                    fprintf(stderr, "%s", str);  /* str includes \n */
            }
        }
    }

  FINISHED:
    if (!found_help)
        fprintf(stderr, "No help found for program <%s>\n", program_name);
    else if (found_help != 2)
        fprintf(stderr, "Warning: no matching tag found\n");
    else {
        fprintf(stderr, "AUTHOR\n%s%s\n\n", prefix, Gvars.X3DNA_VER);
        fprintf(stderr, "Please post questions/comments on the 3DNA Forum:"
                " http://forum.x3dna.org/\n");
        fprintf(stderr, "Please check 'http://x3dna.org/citations' on how to cite"
                " 3DNA --- THANKS!\n");
    }

    print_sep(stderr, '=', nchar);

    close_file(fp);
}

/* Delete all H-atoms and write the rest coordinates in PDB format */
void delH_pdbfile(char *inpfile, char *outfile)
{
    char *ChainID, *p, **AtomName, **ResName, **Miscs;
    double **xyz;
    long i, k = 0, num, *ResSeq, *idx;
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

    idx = lvector(1, num);
    atom_idx(num, AtomName, Miscs, idx);  /* atom index to mark H */

    fp = open_file(outfile, "w");
    print_pdb_title(inpfile, "*", fp);
    for (i = 1; i <= num; i++) {
        if (idx[i] == 3)  /* H-atom */
            continue;
        p = Miscs[i];
        fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n",
                (p[0] == 'A') ? "ATOM  " : "HETATM", ++k, AtomName[i],
                p[1], ResName[i], ChainID[i], ResSeq[i], p[2], xyz[i][1],
                xyz[i][2], xyz[i][3], p + 3);
    }
    fprintf(fp, "END\n");
    close_file(fp);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lvector(idx, 1, num);
}

void contact_msg(long prt_msg)
{
    if (prt_msg) {
        fprintf(stderr, "%s\n", Gvars.X3DNA_VER);
        help3dna("contact_info");
    }
    exit(1);
}

/* partial command-line option match from the beginning */
long str_pmatch(char *str, char *sstr)
{
    return !strncmp(str, sstr, strlen(sstr));
}

/* case-insensitive command-line option match from the beginning */
long case_str_pmatch(char *str, char *sstr)
{
    return !case_strncmp(str, sstr, strlen(sstr));
}

static long isEMPTY_or_isNAN_or_isINF(char *str)
{
    lowerstr(str);  /* to lower case for consistency: NAN/nan etc */

    /* empty field/missing value etc: \t\t; NaN; Inf/-Inf */
    if (*str == '\0' || strstr(str, "nan") || strstr(str, "inf"))
        return TRUE;
    else
        return FALSE;
}

/* check if string "str" contains only valid numerical values */
long is_numeric(char *str)
{
    char *endp;
    double d;

    errno = 0;
    d = strtod(str, &endp);
    UNUSED_PARAMETER(d);

    if (*endp != '\0' || errno == ERANGE)
        return 0L;
    else
        return 1L;
}

/* convert string "str" to a double value, with error checking */
double cvt2double(char *str)
{
    char *endp, *p = trim(str);
    double d;

    if (isEMPTY_or_isNAN_or_isINF(p))
        return XBIG;  /* as an impossible value */

    errno = 0;
    d = strtod(p, &endp);
    if (*endp != '\0' || errno == ERANGE)
        return XBIG;  /* for NA etc */
    else
        return d;
}

/* convert string "str" to a long value, with error checking */
long cvt2long(char *str)
{
    char *endp, *p = trim(str);
    long d;

    if (isEMPTY_or_isNAN_or_isINF(p))
        return LONG_MAX;  /* as an impossible value */

    errno = 0;
    d = strtol(p, &endp, 10);
    if (*endp != '\0' || errno == ERANGE)
        return LONG_MAX;  /* for NA etc */
    else
        return d;
}

long equalsign_pos(char *str)
{
    char *ptr;

    ptr = strchr(str, '=');  /* check from the beginning */
    if (ptr == NULL)
        fatal("wrong format for var=value pair [%s]\n", str);
    return ptr - str + 1;
}

/* extract a long numerical value from command line */
long get_lvalue(char *str, long vmin, long vmax)
{
    long npos, val;

    npos = equalsign_pos(str);

    val = cvt2long(str + npos);
    if (val == LONG_MAX)
        fatal("wrong option [%s]: not a valid integer value\n", str);

    if (!lval_in_range(val, vmin, vmax))
        fatal("invalid option [%s]: value %ld out of range [%ld %ld]\n", str, val, vmin, vmax);

    return val;
}

/* extract a double numerical value from command line */
double get_dvalue(char *str, double vmin, double vmax)
{
    long npos;
    double val;

    npos = equalsign_pos(str);

    val = cvt2double(str + npos);
    if (val > XBIG_CUTOFF)
        fatal("wrong option [%s]: not a valid numerical value\n", str);

    if (!dval_in_range(val, vmin, vmax))
        fatal("invalid option [%s]: value %f out of range [%f %f]\n", str, val, vmin, vmax);

    return val;
}

/* extract a string from command line, substitute ~ to $HOME */
void get_strvalue(char *str, char *dst, long expand_tilde)
{
    char *p;
    long npos;

    if (strlen(str) > BUF512)
        fatal("command line option too long: %.36s...\n", str);

    npos = equalsign_pos(str);
    if (expand_tilde && str[npos] == '~') {
        p = getenv("HOME");
        if (p == NULL) {
            fprintf(stderr, "no environment variable HOME defined!\n");
            strcpy(dst, str + npos);
        } else
            sprintf(dst, "%s%s", p, str + npos + 1);
    } else
        strcpy(dst, str + npos);
}

/* reverse string 'str' in place */
void reverse_string(char *str)
{
    long i = 0, j = strlen(str) - 1;

    while (i < j) {
        cval_swap(&str[i], &str[j]);
        i++;
        j--;
    }
}

/* convert all occurrences of any character in set1 to 'c2' in string 'str' */
void cvtstr_set1toc2(char *str, char *set1, char c2)
{
    char *p = str;

    while (*p) {
        if (strchr(set1, *p))
            *p = c2;
        p++;
    }
}

/* convert all occurrences of character 'c1' to 'c2' in string 'str' */
void cvtstr_c1toc2(char *str, char c1, char c2)
{
    char *p = str;

    while (*p) {
        if (*p == c1)
            *p = c2;
        p++;
    }
}

/* assuming z1 and z2 are already normalized */
double z1_z2_angle_in_0_to_90(double *z1, double *z2)
{
    double dircos = dot(z1, z2);

    return 90.0 - fabs(dot2ang(dircos) - 90.0);  /* absolute value */
}

void do_nothing(void)
{
    return;
}

void skip_lines(long num, FILE * fp)
{
    char *p0;
    long i;

    for (i = 1; i <= num; i++) {
        p0 = my_getline(fp);
        if (p0 == NULL)
            fatal("no <%ld> lines to skip!\n", num);
        free(p0);
    }
}

void check_havefile(char *filename, char *msg)
{
    if (filename == NULL || is_empty_string(filename))
        fatal("File name not specified [%s]\n", msg);

    if (!exist_file(filename))
        fatal("Specified file <%s> does not exist [%s]\n", filename, msg);
}
