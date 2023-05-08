#include "x3dna.h"

/* standard error handler */
void nrerror(char *error_text)
{
    fprintf(stderr, "%s\n", error_text);
    exit(1);
}

void vector_boundary_check(long nl, long nh, char *fun_name)
{
    return;  /* not used */

    if (nl > nh)
        fprintf(stderr, "boundary for %s: low = %ld > high = %ld\n", fun_name, nl, nh);
}

void matrix_boundary_check(long nrl, long nrh, long ncl, long nch, char *fun_name)
{
    return;  /* not used */

    if (nrl > nrh || ncl > nch)
        fprintf(stderr, "boundary for %s: [%ld to %ld; %ld to %ld]\n",
                fun_name, nrl, nrh, ncl, nch);
}

/* ------------------------------------------------------------------ */
/* allocate a char vector with subscript range v[nl..nh] */
char *cvector(long nl, long nh)
{
    char *v;

    vector_boundary_check(nl, nh, "cvector()");
    v = cvector_nr(nl, nh);
    init_cvector(v, nl, nh, '\0');

    return v;
}

/* allocate a char vector with subscript range v[nl..nh] */
char *cvector_nr(long nl, long nh)
{
    char *v;

    v = (char *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(char)));
    if (!v)
        nrerror("allocation failure in cvector()");

    return v - nl + NR_END;
}

/* ------------------------------------------------------------------ */
/* allocate a double vector with subscript range v[nl..nh] */
double *dvector(long nl, long nh)
{
    double *v;

    vector_boundary_check(nl, nh, "dvector()");
    v = dvector_nr(nl, nh);
    init_dvector(v, nl, nh, 0.0);

    return v;
}

/* allocate a double vector with subscript range v[nl..nh] */
double *dvector_nr(long nl, long nh)
{
    double *v;

    v = (double *)
        malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v)
        nrerror("allocation failure in dvector()");

    return v - nl + NR_END;
}

/* ------------------------------------------------------------------ */
/* allocate a long vector with subscript range v[nl..nh] */
long *lvector(long nl, long nh)
{
    long *v;

    vector_boundary_check(nl, nh, "lvector()");
    v = lvector_nr(nl, nh);
    init_lvector(v, nl, nh, 0);

    return v;
}

/* allocate a long vector with subscript range v[nl..nh] */
long *lvector_nr(long nl, long nh)
{
    long *v;

    v = (long *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(long)));
    if (!v)
        nrerror("allocation failure in lvector()");

    return v - nl + NR_END;
}

/* ------------------------------------------------------------------ */
/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
char **cmatrix(long nrl, long nrh, long ncl, long nch)
{
    char **m;

    matrix_boundary_check(nrl, nrh, ncl, nch, "cmatrix()");
    m = cmatrix_nr(nrl, nrh, ncl, nch);
    init_cmatrix(m, nrl, nrh, ncl, nch, '\0');

    return m;
}

/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
char **cmatrix_nr(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    char **m;

    /* allocate pointers to rows */
    m = (char **) malloc((size_t) ((nrow + NR_END) * sizeof(char *)));
    if (!m)
        nrerror("allocation failure 1 in cmatrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (char *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(char)));
    if (!m[nrl])
        nrerror("allocation failure 2 in cmatrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

/* ------------------------------------------------------------------ */
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
    double **m;

    matrix_boundary_check(nrl, nrh, ncl, nch, "dmatrix()");
    m = dmatrix_nr(nrl, nrh, ncl, nch);
    init_dmatrix(m, nrl, nrh, ncl, nch, 0.0);

    return m;
}

/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix_nr(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;

    /* allocate pointers to rows */
    m = (double **) malloc((size_t) ((nrow + NR_END) * sizeof(double *)));
    if (!m)
        nrerror("allocation failure 1 in dmatrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double *)
        malloc((size_t) ((nrow * ncol + NR_END) * sizeof(double)));
    if (!m[nrl])
        nrerror("allocation failure 2 in dmatrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

/* ------------------------------------------------------------------ */
/* allocate a long matrix with subscript range m[nrl..nrh][ncl..nch] */
long **lmatrix(long nrl, long nrh, long ncl, long nch)
{
    long **m;

    matrix_boundary_check(nrl, nrh, ncl, nch, "lmatrix()");
    m = lmatrix_nr(nrl, nrh, ncl, nch);
    init_lmatrix(m, nrl, nrh, ncl, nch, 0);

    return m;
}

/* allocate a long matrix with subscript range m[nrl..nrh][ncl..nch] */
long **lmatrix_nr(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    long **m;

    /* allocate pointers to rows */
    m = (long **) malloc((size_t) ((nrow + NR_END) * sizeof(long *)));
    if (!m)
        nrerror("allocation failure 1 in lmatrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (long *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(long)));
    if (!m[nrl])
        nrerror("allocation failure 2 in lmatrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

/* ------------------------------------------------------------------ */
/* free a char vector allocated with cvector() */
void free_cvector(char *v, long nl, long nh)
{
    free((FREE_ARG) (v + nl - NR_END));
    UNUSED_PARAMETER(nh);  /* strict compiler options */
}

/* free a double vector allocated with dvector() */
void free_dvector(double *v, long nl, long nh)
{
    free((FREE_ARG) (v + nl - NR_END));
    UNUSED_PARAMETER(nh);  /* strict compiler options */
}

/* free a long vector allocated with lvector() */
void free_lvector(long *v, long nl, long nh)
{
    free((FREE_ARG) (v + nl - NR_END));
    UNUSED_PARAMETER(nh);  /* strict compiler options */
}

/* free a char matrix allocated by cmatrix() */
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    UNUSED_PARAMETER(nrh);
    UNUSED_PARAMETER(nch);  /* strict compiler options */
}

/* free a double matrix allocated by dmatrix() */
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    UNUSED_PARAMETER(nrh);
    UNUSED_PARAMETER(nch);  /* strict compiler options */
}

/* free a long matrix allocated by lmatrix() */
void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    UNUSED_PARAMETER(nrh);
    UNUSED_PARAMETER(nch);  /* strict compiler options */
}

/* ------------------------------------------------------------------ */
double dval_sqr(double dval)
{
    if (dval == 0.0)
        return 0.0;
    else
        return dval * dval;
}

void dval_swap(double *pa, double *pb)
{
    double temp;

    temp = *pa;
    *pa = *pb;
    *pb = temp;
}

void lval_swap(long *pa, long *pb)
{
    long temp;

    temp = *pa;
    *pa = *pb;
    *pb = temp;
}

void cval_swap(char *pa, char *pb)
{
    int c;

    c = *pa;
    *pa = *pb;
    *pb = c;
}

double dval_max(double a, double b)
{
    return (a > b) ? a : b;
}

double dval_min(double a, double b)
{
    return (a < b) ? a : b;
}

long lval_max(long a, long b)
{
    return (a > b) ? a : b;
}

long lval_min(long a, long b)
{
    return (a < b) ? a : b;
}

/* absolute difference between two doubles */
double abs_dval_diff(double a, double b)
{
    return fabs(a - b);
}

/* check if lval exists in s[ib] .. s[ie] */
long lval_in_set(long lval, long ib, long ie, long *s)
{
    long i;

    for (i = ib; i <= ie; i++)
        if (lval == s[i])
            return TRUE;

    return FALSE;
}

/* check if "dval" in within "dlow" and "dhigh" */
long dval_in_range(double dval, double dlow, double dhigh)
{
    return (dval >= dlow && dval <= dhigh) ? 1L : 0L;
}

/* check if "lval" in within "llow" and "lhigh" */
long lval_in_range(long lval, long llow, long lhigh)
{
    return (lval >= llow && lval <= lhigh) ? 1L : 0L;
}

void max_dmatrix(double **d, long nr, long nc, double *maxdm)
{
    long i, j;

    for (i = 1; i <= nc; i++) {
        maxdm[i] = -XBIG;
        for (j = 1; j <= nr; j++)
            maxdm[i] = dval_max(maxdm[i], d[j][i]);
    }
}

void min_dmatrix(double **d, long nr, long nc, double *mindm)
{
    long i, j;

    for (i = 1; i <= nc; i++) {
        mindm[i] = XBIG;
        for (j = 1; j <= nr; j++)
            mindm[i] = dval_min(mindm[i], d[j][i]);
    }
}

void ave_dmatrix(double **d, long nr, long nc, double *avedm)
{
    long i, j;

    if (!nr)
        nrerror("divided by zero in <ave_dmatrix>");

    for (i = 1; i <= nc; i++) {
        avedm[i] = 0.0;
        for (j = 1; j <= nr; j++)
            avedm[i] += d[j][i];
        avedm[i] /= nr;
    }
}

void std_dmatrix(double **d, long nr, long nc, double *stddm)
{
    double dsum, temp;
    double *aved;
    long i, j;

    if (nr < 2) {
        fprintf(stderr, "less than two items for <std_dmatrix>\n");
        init_dvector(stddm, 1, nc, 0.0);  /* column-wise */
        return;
    }

    aved = dvector(1, nc);
    ave_dmatrix(d, nr, nc, aved);

    for (i = 1; i <= nc; i++) {
        dsum = 0.0;
        for (j = 1; j <= nr; j++) {
            temp = d[j][i] - aved[i];
            dsum += temp * temp;
        }
        stddm[i] = sqrt(dsum / (nr - 1));
    }

    free_dvector(aved, 1, nc);
}

double max_dvector(double *d, long nl, long nh)
{
    double maxdv = -XBIG;
    long i;

    for (i = nl; i <= nh; i++)
        maxdv = dval_max(maxdv, d[i]);

    return maxdv;
}

double min_dvector(double *d, long nl, long nh)
{
    double mindv = XBIG;
    long i;

    for (i = nl; i <= nh; i++)
        mindv = dval_min(mindv, d[i]);

    return mindv;
}

double ave_dvector(double *d, long n)
{
    double dsum = 0.0;
    long i;

    if (n <= 0) {
        fprintf(stderr, "zero (or negative) items [%ld] for <ave_dvector>\n", n);
        return 0.0;
    }

    for (i = 1; i <= n; i++)
        dsum += d[i];

    return dsum / n;
}

double std_dvector(double *d, long n)
{
    double aved, dsum = 0.0, temp;
    long i;

    if (n < 2) {
        fprintf(stderr, "less than 2 items for <std_dvector>\n");
        return 0.0;
    }

    aved = ave_dvector(d, n);
    for (i = 1; i <= n; i++) {
        temp = d[i] - aved;
        dsum += temp * temp;
    }

    return sqrt(dsum / (n - 1));
}

/* initialize character-matrix dmtx with init_val */
void init_cmatrix(char **cmtx, long nrl, long nrh, long ncl, long nch, char init_val)
{
    long i;

    for (i = nrl; i <= nrh; i++)
        init_cvector(cmtx[i], ncl, nch, init_val);
}

/* initialize double-matrix dmtx with init_val */
void init_dmatrix(double **dmtx, long nrl, long nrh, long ncl, long nch, double init_val)
{
    long i;

    for (i = nrl; i <= nrh; i++)
        init_dvector(dmtx[i], ncl, nch, init_val);
}

/* initialize long-matrix lmtx with init_val */
void init_lmatrix(long **lmtx, long nrl, long nrh, long ncl, long nch, long init_val)
{
    long i;

    for (i = nrl; i <= nrh; i++)
        init_lvector(lmtx[i], ncl, nch, init_val);
}

/* initialize a character vector cvec[ib..ie] with init_val, as a string */
void init_cvector(char *cvec, long ib, long ie, char init_val)
{
    long i;

    for (i = ib; i < ie; i++)
        cvec[i] = init_val;
    cvec[ie] = '\0';
}

/* initialize a character vector cvec[ib..ie] with init_val, as a vector */
void init_cvector_all(char *cvec, long ib, long ie, char init_val)
{
    long i;

    for (i = ib; i <= ie; i++)
        cvec[i] = init_val;
}

/* initialize a double vector dvec[ib..ie] with init_val */
void init_dvector(double *dvec, long ib, long ie, double init_val)
{
    long i;

    for (i = ib; i <= ie; i++)
        dvec[i] = init_val;
}

/* initialize a long vector lvec[ib..ie] with init_val */
void init_lvector(long *lvec, long ib, long ie, long init_val)
{
    long i;

    for (i = ib; i <= ie; i++)
        lvec[i] = init_val;
}

/* copy a double vector [nl .. nh], from s to d */
void copy_dvector(double *d, double *s, long nl, long nh)
{
    long i;

    for (i = nl; i <= nh; i++)
        d[i] = s[i];
}

/* copy a long vector [nl .. nh], from s to d */
void copy_lvector(long *d, long *s, long nl, long nh)
{
    long i;

    for (i = nl; i <= nh; i++)
        d[i] = s[i];
}

int dval_compare(const void *v1, const void *v2)
{
    const double *p1, *p2;

    p1 = (const double *) v1;
    p2 = (const double *) v2;

    return (*p1 > *p2) ? 1 : (*p1 < *p2) ? -1 : 0;
}

int lval_compare(const void *v1, const void *v2)
{
    const long *p1, *p2;

    p1 = (const long *) v1;
    p2 = (const long *) v2;

    return (*p1 > *p2) ? 1 : (*p1 < *p2) ? -1 : 0;
}

int cstr_compare(const void *v1, const void *v2)
{
    const char **p1, **p2;

    p1 = (const char **) v1;
    p2 = (const char **) v2;

    return strcmp(*p1, *p2);
}

/* negate each value of vector xyz1 */
void negate_xyz(double *xyz1)
{
    long i;

    for (i = 1; i <= 3; i++)
        xyz1[i] = -xyz1[i];
}

/* distance between vectors xyz1 and xyz2 */
double p1p2_dist(double *xyz1, double *xyz2)
{
    double dxyz[4];

    ddxyz(xyz1, xyz2, dxyz);
    return veclen(dxyz);
}

void p1p2_ave(double *xyz1, double *xyz2, double *ave)
{
    long i;

    for (i = 1; i <= 3; i++)
        ave[i] = 0.5 * (xyz1[i] + xyz2[i]);
}

/* check if distance between xyz1 & xyz2 is within (dlow, dhigh) */
long within_limits(double *xyz1, double *xyz2, double dlow, double dhigh)
{
    double d, dxyz[4];
    long i, ltok = 0;

    for (i = 1; i <= 3; i++)
        if (fabs(dxyz[i] = xyz1[i] - xyz2[i]) > dhigh)
            break;
    if (i > 3) {
        d = veclen(dxyz);
        if (dval_in_range(d, dlow, dhigh))
            ltok = 1;
    }
    return ltok;
}

/* get the sum vector of xyz1 and xyz2 */
void sumxyz(double *xyz1, double *xyz2, double *sxyz)
{
    long i;

    for (i = 1; i <= 3; i++)
        sxyz[i] = xyz1[i] + xyz2[i];
}

/* get the mean vector of xyz1 and xyz2 */
void avexyz(double *xyz1, double *xyz2, double *mxyz)
{
    long i;

    for (i = 1; i <= 3; i++)
        mxyz[i] = 0.5 * (xyz1[i] + xyz2[i]);
}

/* the difference vector of xyz2 and xyz1 */
void ddxyz(double *xyz1, double *xyz2, double *dxyz)
{
    long i;

    for (i = 1; i <= 3; i++)
        dxyz[i] = xyz2[i] - xyz1[i];
}

/* make a copy of xyz1 to xyz2 */
void cpxyz(double *xyz1, double *xyz2)
{
    long i;

    for (i = 1; i <= 3; i++)
        xyz2[i] = xyz1[i];
}

/* get orthogonal component of va w.r.t. vref [1-by-3] */
void vec_orth(double *va, double *vref)
{
    double d;
    long i;

    vec_norm(vref);
    d = dot(va, vref);
    for (i = 1; i <= 3; i++)
        va[i] -= d * vref[i];
    vec_norm(va);
}

/* dot product between two 1-by-3 vectors */
double dot(double *va, double *vb)
{
    double dsum = 0.0;
    long i;

    for (i = 1; i <= 3; i++)
        dsum += va[i] * vb[i];

    return dsum;
}

/* cross product between two 1-by-3 vectors */
void cross(double *va, double *vb, double *vc)
{
    vc[1] = va[2] * vb[3] - va[3] * vb[2];
    vc[2] = va[3] * vb[1] - va[1] * vb[3];
    vc[3] = va[1] * vb[2] - va[2] * vb[1];
}

long sign_control(double *va, double *vb, double *vref)
{
    double dtmp[4];

    cross(va, vb, dtmp);
    if (dot(dtmp, vref) < 0.0)
        return -1;
    else
        return 1;
}

/* length (magnitude) of a 1-by-3 vector */
double veclen(double *va)
{
    return sqrt(dot(va, va));
}

/* normalize a 1-by-3 vector, i.e. with unit length */
void vec_norm(double *va)
{
    double vlen;
    long i;

    vlen = veclen(va);
    if (vlen > XEPS)
        for (i = 1; i <= 3; i++)
            va[i] /= vlen;
}

/* from dot product of 2 NORMAL vectors to angle: positive */
double dot2ang(double dotval)
{
    double ang_deg;
    if (dotval >= 1.0)
        ang_deg = 0.0;
    else if (dotval <= -1.0)
        ang_deg = 180.0;
    else
        ang_deg = rad2deg(acos(dotval));
    return ang_deg;
}

/* angle magnitude in degrees between two 1-by-3 vectors */
double magang(double *va, double *vb)
{
    double ang_deg;
    if (veclen(va) < XEPS || veclen(vb) < XEPS)
        ang_deg = 0.0;
    else {
        vec_norm(va);
        vec_norm(vb);
        ang_deg = dot2ang(dot(va, vb));
    }
    return ang_deg;
}

double rad2deg(double ang)
{
    return ang * 180.0 / PI;
}

double deg2rad(double ang)
{
    return ang * PI / 180.0;
}

void copy_dmatrix(double **a, long nr, long nc, double **o)
{
    long i, j;

    for (i = 1; i <= nr; i++)
        for (j = 1; j <= nc; j++)
            o[i][j] = a[i][j];
}

void copy_lmatrix(long **a, long nr, long nc, long **o)
{
    long i, j;

    for (i = 1; i <= nr; i++)
        for (j = 1; j <= nc; j++)
            o[i][j] = a[i][j];
}

void multi_matrix(double **a, long nra, long nca, double **b, long nrb, long ncb, double **o)
{
    long i, j, k;

    if (nca != nrb)
        nrerror("matrices a and b do not conform");

    for (i = 1; i <= nra; i++) {
        for (j = 1; j <= ncb; j++) {
            o[i][j] = 0.0;
            for (k = 1; k <= nca; k++)
                o[i][j] += a[i][k] * b[k][j];
        }
    }
}

/* vector-matrix multiplication */
void multi_vec_matrix(double *a, long n, double **b, long nr, long nc, double *o)
{
    long i, j;

    if (n != nr)
        nrerror("vector and matrix do not conform");

    for (i = 1; i <= nc; i++) {
        o[i] = 0.0;
        for (j = 1; j <= n; j++)
            o[i] += a[j] * b[j][i];
    }
}

/* vector - transpose-of-matrix multiplication */
void multi_vec_Tmatrix(double *a, long n, double **b, long nr, long nc, double *o)
{
    double **tb;  /* transpose of b */

    tb = dmatrix(1, nc, 1, nr);

    transpose_matrix(b, nr, nc, tb);
    multi_vec_matrix(a, n, tb, nc, nr, o);

    free_dmatrix(tb, 1, nc, 1, nr);
}

void transpose_matrix(double **a, long nr, long nc, double **o)
{
    long i, j;

    for (i = 1; i <= nc; i++)
        for (j = 1; j <= nr; j++)
            o[i][j] = a[j][i];
}

void identity_matrix(double **d, long n)
{
    long i;

    for (i = 1; i <= n; i++) {
        init_dvector(d[i], 1, n, 0.0);
        d[i][i] = 1.0;
    }
}
