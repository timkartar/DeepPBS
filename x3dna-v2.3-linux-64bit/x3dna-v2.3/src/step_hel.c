#include "x3dna.h"

typedef struct {
    long hel2step;
} struct_args;

static void read_six_pars(double *pars)
{
    char str[BUF512];
    long i, ik;

    while (1)
        if (fgets(str, sizeof str, stdin) != NULL) {
            for (i = 0; i < (long) strlen(str); i++)
                if (str[i] == ',')
                    str[i] = ' ';  /* change comma to space */
            ik = sscanf(str, "%lf %lf %lf %lf %lf %lf", &pars[1], &pars[2],
                        &pars[3], &pars[4], &pars[5], &pars[6]);
            if (ik == 6)
                return;
            else
                fprintf(stderr, "please input six parameters\n");
        } else
            fatal("error reading six parameters\n");
}

static void step_hel_usage(void)
{
    help3dna_usage("step_hel");
}

static void step_hel_cmdline(int argc, char *argv[], struct_args * args)
{
    long i;

    args->hel2step = FALSE;

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        if (str_pmatch(argv[i], "-h"))
            step_hel_usage();
        else if (str_pmatch(argv[1], "-s"))
            args->hel2step = TRUE;
        else {
            fprintf(stderr, "\tunrecognized option [%s]\n", argv[i]);
            step_hel_usage();
        }
    }

    if (argc > i) {
        fprintf(stderr, "\tWrong option [%s]\n", argv[i]);
        step_hel_usage();
    }
}

/* interchange between helical parameters and step parameters */
int main(int argc, char *argv[])
{
    char *format = "%10.4f";
    double inp_par[7], out_par[7], chk_par[7];
    double o1[4] = { EMPTY_NUMBER, 0.0, 0.0, 0.0 }, o2[4], mpos[4], hmpos[4], cmpos[4];
    double **r1, **r2, **mst, **hmst, **cmst;
    long i;
    struct_args args;

    set_my_globals(argv[0]);

    step_hel_cmdline(argc, argv, &args);

    if (args.hel2step)
        fprintf(stderr, "\nSix helical parameters in the order of: \n"
                "X-disp  Y-disp  h-Rise  Incl.  Tip  h-Twist\n");
    else
        fprintf(stderr, "\nSix step parameters in the order of: \n"
                "Shift  Slide  Rise  Tilt  Roll  Twist\n");

    read_six_pars(inp_par);

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);
    hmst = dmatrix(1, 3, 1, 3);
    cmst = dmatrix(1, 3, 1, 3);

    identity_matrix(r1, 3);  /* identity-matrix */
    if (args.hel2step) {  /* helical to step */
        xhelfunc(inp_par, r2, hmst, o2, hmpos);
        bpstep_par(r1, o1, r2, o2, out_par, mst, mpos);
        fprintf(stdout, "\nYour six input helical parameters: \n"
                "    X-disp    Y-disp    h-Rise     Incl.       Tip   h-Twist\n");
        for (i = 1; i <= 6; i++)
            fprintf(stdout, format, inp_par[i]);
        fprintf(stdout, "\n");
        fprintf(stdout, "\nYour six output step parameters: \n"
                "     Shift     Slide      Rise      Tilt      Roll     Twist\n");
        for (i = 1; i <= 6; i++)
            fprintf(stdout, format, out_par[i]);
        fprintf(stdout, "\n");

    } else {  /* step to helical */
        xbpfunc(inp_par, r2, mst, o2, mpos);
        /* re-calculate the step parameters for verification */
        bpstep_par(r1, o1, r2, o2, chk_par, cmst, cmpos);
        helical_par(r1, o1, r2, o2, out_par, hmst, hmpos);
        fprintf(stdout, "\nYour six input step parameters:\n"
                "     Shift     Slide     Rise      Tilt     Roll       Twist\n");
        for (i = 1; i <= 6; i++)
            fprintf(stdout, format, inp_par[i]);
        fprintf(stdout, "\nRecalculated parameters (***for verification***):\n");
        for (i = 1; i <= 6; i++)
            fprintf(stdout, format, chk_par[i]);
        fprintf(stdout, "\n");
        fprintf(stdout, "\nYour six output helical parameters:\n"
                "    X-disp    Y-disp    h-Rise     Incl.     Tip     h-Twist\n");
        for (i = 1; i <= 6; i++)
            fprintf(stdout, format, out_par[i]);
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\nPosition of the 2nd base-pair: %9.4f%9.4f%9.4f\n\n", o2[1], o2[2], o2[3]);
    fprintf(stdout, "Direction cosines of x-axis:   %9.4f%9.4f%9.4f\n",
            r2[1][1], r2[2][1], r2[3][1]);
    fprintf(stdout, "Direction cosines of y-axis:   %9.4f%9.4f%9.4f\n",
            r2[1][2], r2[2][2], r2[3][2]);
    fprintf(stdout, "Direction cosines of z-axis:   %9.4f%9.4f%9.4f\n",
            r2[1][3], r2[2][3], r2[3][3]);

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(hmst, 1, 3, 1, 3);

    clear_my_globals();

    return 0;
}
