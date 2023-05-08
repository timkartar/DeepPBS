#include "x3dna.h"

typedef struct {
    char inppdb[BUF512];
    char outpdb[BUF512];
} struct_args;

static void o1p_o2p_usage(void)
{
    help3dna_usage("o1p_o2p");
}

static void o1p_o2p_cmdline(int argc, char *argv[], struct_args * args)
{
    long i;

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;

        if (check_global_options(argv[i]))
            continue;

        o1p_o2p_usage();  /* no other options allowed */
    }

    if (argc == i + 2) {
        strcpy(args->inppdb, argv[i]);
        strcpy(args->outpdb, argv[i + 1]);
    } else
        o1p_o2p_usage();
}

/* there is a confusion about O1P/O2P labeling in nucleic acid
   structures which affects least-fitting RMS value between two
   structures. This utility checks if O1P/O2P in a PDB file are
   properly labeled as in the NDB. A new clean PDB file is written. */
int main(int argc, char *argv[])
{
    char idmsg[BUF512], *ChainID, **AtomName, **ResName, **Miscs;
    double d, O1P_O2P[4], O2P_O5[4], O_P[4], temp[4], **xyz;
    long i, ib, ie, inum = 0, j, num, num_residue, O1P, O2P, O5, P;
    long *ResSeq, **seidx;
    FILE *fp;
    struct_args args;

    set_my_globals(argv[0]);

    o1p_o2p_cmdline(argc, argv, &args);

    /* read in PDB file */
    num = number_of_atoms(args.inppdb, 1, "*");
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(args.inppdb, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");

    /* get the numbering information of each residue */
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

    /* check if O1P/O2P labeling are okay */
    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];

        if (num_strmatch(" P  ", AtomName, ib, ie)) {  /* with P atom */
            sprintf(idmsg, "residue %3s %4ld%c on chain %c [#%ld]",
                    ResName[ib], ResSeq[ib], Miscs[ib][2], ChainID[ib], i);
            O1P = find_1st_atom(" O1P", AtomName, ib, ie, idmsg);
            O2P = find_1st_atom(" O2P", AtomName, ib, ie, idmsg);
            O5 = find_1st_atom(" O5'", AtomName, ib, ie, idmsg);
            P = find_1st_atom(" P  ", AtomName, ib, ie, idmsg);

            if (O1P && O2P && O5 && P) {
                ddxyz(xyz[O1P], xyz[O2P], O1P_O2P);
                ddxyz(xyz[O2P], xyz[O5], O2P_O5);
                for (j = 1; j <= 3; j++) {
                    d = (xyz[O1P][j] + xyz[O2P][j] + xyz[O5][j]) / 3.0;
                    O_P[j] = xyz[P][j] - d;
                }
                cross(O1P_O2P, O2P_O5, temp);
                if (dot(temp, O_P) < 0) {
                    fprintf(stderr, "*PO4 group at %s does NOT conform to convention*\n", idmsg);
                    for (j = 1; j <= 3; j++)
                        dval_swap(&xyz[O1P][j], &xyz[O2P][j]);
                }
            }
        }
    }
    fp = open_file(args.outpdb, "w");
    print_pdb_title(args.inppdb, "*", fp);
    pdb_record(1, num, &inum, 0, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, fp);
    fprintf(fp, "END\n");
    close_file(fp);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);

    clear_my_globals();

    return 0;
}
