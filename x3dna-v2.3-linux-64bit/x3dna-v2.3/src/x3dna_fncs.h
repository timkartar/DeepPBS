
/* alc2img.c */

/* ana_fncs.c */
void populate_nt_info(long num_residue, long **seidx, char **ResName, char *ChainID,
                      long *ResSeq, char **Miscs, char *bseq, char **nt_info);
void populate_nt_list(long num_residue, long **seidx, long *RY, char *bseq,
                      char **AtomName, double **xyz, long **nt_list);
long **read_input(char *inpfile, char *pdbfile, char *outfile, long *ds, long *num_bp,
                  long *ip, long *hetatm);
void print_header(long ds, long num_bp, long num, char *pdbfile, FILE * fp);
void output_Borg_P_C1_C4(long num_residue, double **org, double **xyz, long **nt_list,
                         char **nt_info);
void atom_list(long ds, long num_bp, long **pair_num, long **seidx, long *RY,
               char **bp_seq, char **AtomName, char **ResName, char *ChainID,
               long *ResSeq, char **Miscs, long **phos, long **c6_c8, long **sugar, long **chi);
void get_nt_torsion(long num_residue, double **org, double **xyz, long **nt_list,
                    double **nt_torsion);
void get_ss_Zp_Dp(long num_residue, double **org, double **orien, double **xyz,
                  long **nt_list, double **ss_Zp_Dp);
void output_nt_torsion(long num_residue, char **nt_info, long **nt_list,
                       double **nt_torsion, double **ss_Zp_Dp, FILE * fp);
void get_nt_bb_torsion(double **nt_bb_torsion, long num_residue, long **seidx,
                       long *RY, char **AtomName, char **ResName, char *ChainID,
                       long *ResSeq, char **Miscs, double **xyz);
void backbone_torsion(long ds, long num_bp, long **pair_num, char **bp_seq, long **sugar,
                      long **chi, double **xyz, double **nt_bb_torsion, FILE * fp);
void p_c1_dist(long ds, long num_bp, char **bp_seq, long **phos, long **chi,
               double **xyz, long *bphlx, FILE * fp);
void lambda_d3(long num_bp, char **bp_seq, long **chi, long **c6_c8, double **xyz, FILE * fp);
void print_axyz(long num_bp, char **bp_seq, long **aidx, char *aname, double **xyz);
void groove_width(long parallel, long num_bp, char **bp_seq, long **phos,
                  double **xyz, long *bphlx, FILE * fp);
void check_wc_wobble_pair(long *bpid, char *bp, double shear, double stretch, double opening);
void set_chain_nmarkers019_to_symbols(long num, long *nmarkers, char *cmarkers);
void get_bp_3char_symbols(long bp_type, char zdir, char *bp_sym);
void ref_frames(long ds, long num_bp, long **pair_num, char **bp_seq, long **seidx,
                long *RY, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                char **Miscs, double **xyz, FILE * fp, double **orien, double **org,
                long *WC_info, long *str_type, long irna, long **o3p_brk);
void bpstep_par(double **rot1, double *org1, double **rot2, double *org2, double *pars,
                double **mst_orien, double *mst_org);
void helical_par(double **rot1, double *org1, double **rot2, double *org2, double *pars,
                 double **mst_orien, double *mst_org);
void print_par(char **bp_seq, long num_bp, long ich, long ishel, double **param, FILE * fp);
void output_ave_std(long num, double **parcln, int dnum, char *fmt, FILE * fp);
void prt_stepstr(char **step_str, long num_step, long *bphlx, long ishel, double **param,
                 FILE * fp);
void prt_step_par(char **bp_seq, long num_bp, long *bphlx, long ishel, double **param, FILE * fp);
void bz_check(double **r1, double *o1, double **r2, double *o2, long bz,
              long *bz_junction, long *z_step);
void get_mtwist(long nbpm1, long *bphlx, long *WC_info, double **twist_rise,
                double *twist_p, double *twist_n);
void get_parameters(long ds, long num_bp, char **bp_seq, double **orien, double **org,
                    long *WC_info, FILE * fp, double **twist_rise, double *mst_orien,
                    double *mst_org, double *mst_orienH, double *mst_orgH, long *bphlx,
                    long istart, long istep, long bz, long *str_type, long **pair_num,
                    char **nt_info);
void parvec2mtx(double *parvec, long num, double **parmtx);
void print_ss_rebuild_pars(double **pars, long num_bp, char *str, char **bp_seq, FILE * fp);
void print_ds_rebuild_pars(double **bp_par, double **step_par, long num_bp, char *str,
                           char **bp_seq, FILE * fp);
void print_ref(char **bp_seq, long num_item, long ich, double *org, double *orien, FILE * fp);
void write_mst(long ds, long num_bp, long **pair_num, char **bp_seq, double *mst_orien,
               double *mst_org, long **seidx, char **AtomName, char **ResName,
               char *ChainID, long *ResSeq, double **xyz, char **Miscs,
               long **htm_water, double **twist_rise, char *strfile);
void print_xyzP(long parallel, long nbpm1, char **bp_seq, long **phos, double *mst_orien,
                double *mst_org, double **xyz, FILE * fp, char *title_str, double **aveP,
                long p_offset);
void print_PP(long parallel, double **twist_rise, long num_bp, char **bp_seq, long **phos,
              double *mst_orien, double *mst_org, double *mst_orienH, double *mst_orgH,
              double **xyz, long *WC_info, long *bphlx, long abi, long **chi, FILE * fp);
void str_classify(double twist_p, double twist_n, long str_type, long parallel,
                  long num_bp, FILE * fp);
double a_hlxdist(long idx, double **xyz, double *hlx_axis, double *hlx_pos);
void print_radius(char **bp_seq, long nbpm1, long ich, double **p_radius,
                  double **o4_radius, double **c1_radius, long *bphlx, FILE * fp);
void helix_radius(long ds, long num_bp, char **bp_seq, double **orien, double **org,
                  long **phos, long **chi, double **xyz, long *bphlx, FILE * fp);
void print_shlx(char **bp_seq, long nbpm1, long ich, double *shlx_orien,
                double *shlx_org, FILE * fp);
void get_helix_axis(long ds, long num_bp, char **bp_seq, double **orien,
                    double **org, long *bphlx, FILE * fp);
void get_axis(long nvec, long **idx, long num, double **xyz, long nb, long *C1b,
              long *C1e, double *std_rise, double *hrise, double *haxis,
              double *hstart, double *hend);
void print_poc_r3d(double *rave, double *hstart, double *hend);
void global_analysis(long ds, long num_bp, long num, char **bp_seq, long **chi,
                     long **phos, double **xyz, FILE * fp);
void base_overlap(long ds, long num_bp, long num, long num_residue, long **pair_num,
                  long *bRY, char **bp_seq, long **seidx, char **AtomName, double **xyz,
                  long *idx, double **orien, double **org, FILE * fp);
long ratom_xyz(long *ratom_list, long only_ring, double **xyz, double *oave, double **oxyz);
void get_zoave(long istep, long ds, double **orien, double **org, double *oave, double *zave);
void get_bp_zoave(long ia, long ib, double **orien, double **org, double *oave, double *zave);
void ring_oidx(long num, long num_residue, long *RY, long **seidx, char **AtomName,
               double **xyz, long *idx, long **ring_atom);
void get_cntatom(long *ringlist, long **connect, long *idx);
double get_oarea(long r1, long r2, long **ring_atom, double *oave, double *zave,
                 double **xyz, long only_ring);
void verify_oarea(void);
void cehs_base_atoms(char **AtomName, long ib, long ie, long *num_batom, long *batom);
void cehs_bppar(double **rot1, double *org1, double **rot2, double *org2, double *pars,
                double **mst_orien, double *mst_org);
void cehs_pars(long num_bp, long istart, long istep, long **pair_num, char **bp_seq,
               long **seidx, long **c6_c8, long *RY, char **AtomName, char **ResName,
               char *ChainID, long *ResSeq, char **Miscs, double **xyz, double *bp_orien,
               double *bp_org, long bz, FILE * fp);
void schnaap_global(long num_bp, long num, char **bp_seq, long **chi, double **xyz,
                    double *bp_orien, double *bp_org, FILE * fp);
void out_cehs(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org, FILE * fp);
void compdna(double **rot1, double *org1, double **rot2, double *org2, double *pars,
             double **mst_orien, double *mst_org);
void out_compdna(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org, FILE * fp);
void my_curves(double **rot1, double *org1, double **rot2, double *org2, double *pars);
void curves_mbt(long ibp, double **orien, double **org, double **cvr, double *cvo);
void out_curves(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org, FILE * fp);
void freehelix(double **rot1, double *org1, double **rot2, double *org2, double *pars,
               double **mst_orien, double *mst_org);
void out_freehelix(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org,
                   FILE * fp);
void sgl_helix(double **rot1, double **rot2, double *rot_ang, double *rot_hlx);
void ngeom(double **rot1, double *org1, double **rot2, double *org2, double *pars,
           double **mst_orien, double *mst_org);
void out_ngeom(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org, FILE * fp);
void nuparm(double **rot1, double *org1, double **rot2, double *org2, double *pars,
            double **mst_orien, double *mst_org, double *hpars, long get_hpar);
void out_nuparm(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org, FILE * fp);
void rna(double **rot1, double *org1, double *pvt1, double **rot2, double *org2,
         double *pvt2, double *pars, double **mst_orien, double *mst_org);
void pvt_dxdy(double **rot1, double *org1, double *pvt1, double *pars,
              double **mst_orien, double *mst_org);
void out_rna(long ds, long num_bp, char **bp_seq, long *bphlx, double **orien,
             double **org, FILE * fp);
void other_pars(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org);

/* analyze.c */

/* anyhelix.c */

/* app_fncs.c */
long string_contains_only_those_characters(char *str, char *chars_set);
void bpid_wc_str(long bpid, double zdir, char *wc);
int case_strcmp(const char *s1, const char *s2);
int case_strncmp(const char *s1, const char *s2, long n);
char *case_strstr(const char *haystack, const char *needle);
char *case_strchr(const char *s, int c);
void kbd_input(char *msg);
long get_line_number(char *filename, long skips);
long is_empty_string(const char *str);
long is_equal_string(const char *str1, const char *str2);
long is_equal_case_string(const char *str1, const char *str2);
void null_line_comment(char *str);
long is_comment_line(char *line);
long is_empty_line(char *line);
long is_skip_line(char *line);
void bname_ext(char *src, char *ext, char *dst);
double get_point2line_perp_distance(double *pnt, double *line_p1, double *line_p2);
void get_tag_string_pair(char *prefix, char *tag, char *btag, char *etag);
void get_xml_tag(FILE * fpxml, char *prefix, char *line, char *connector, char *tag, char *tag_str);
void print_xml_tag(FILE * fpxml, char *prefix, char *line, char *tag, char *otag, FILE * fp);
void get_xml_tag_long(FILE * fpxml, char *prefix, char *line, char *tag, long *lval);
void get_xml_tag_double(FILE * fpxml, char *prefix, char *line, char *tag, double *dval);
long tag_match(char *prefix, char *line, char *tag);
void extract_attribute(char *line, char *attr, char *attr_val, long to_lower);
void print_ptable();
long set_switch_default_true(char *option);
void check_required_file(char *filename, char *invalid_str, char *msg);
void define_frame_by_3atoms(long num, double **xyz, double **refmat);
long get_xmlArgNumber(char *str, char *pstr);
void base_frame(long num_residue, char *bseq, long **seidx, long *res_type,
                char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                char **Miscs, double **xyz, char *BDIR, double **orien, double **org);
void peptide_info(long num_residue, long **seidx, char **AtomName, char **ResName,
                  char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                  long *res_type, long *cidx, long *mchain);
void base_blks(long num_residue, long *res_type, double **orien, double **org,
               char *bseq, char *BDIR, char *alcfile);
void set_default_misc_pars(miscPars * misc_pars);
void lsplane_xyz(double **xyz, long num_plane, long *atom_plane, double **nxyz, double *z);
long read_PairInfo(char *inpfile, long **pair_info);
long is_linked_by_gap(long i, long j, double **o3_p);
long is_linked(long i, long j, double **o3_p);
double distance_ab(double **o3_p, long ia, long ib, long ipa, long ipb);
void write_rotmat(double **rotmat);
void read_rotmat(char *rotfile, double **rotmat);
void reverse_y_z_columns(double **R);
long get_num_nt(long num_residue, long *RY);
void peptide_frame(long num_residue, char *BDIR, long *res_type, long *mchain,
                   double **xyz, double **orien, double **org);
void peptide_blks(long num_residue, char *BDIR, long *cidx, double **orien, double **org,
                  char *alcfile);
void set_my_globals(char *pgname);
void clear_my_globals(void);
long check_global_options(char *option);
void get_AA_frames(long num_residue, long **seidx, long *res_type, char **AtomName,
                   char **ResName, char *ChainID, long *ResSeq, char **Miscs,
                   double **xyz, char *BDIR, double **orien, double **org);
void verify_Cb_coordinates(long nmatch, double **ePxyz, double **fitted_xyz);
void residue_chain_resnum(char chain_id, long res_seq, char *misc, char *idmsg);
void residue_strid(char chain_id, long res_seq, char *misc, char *rname, char *idmsg);
void convert_resNameSpace(char *resName, char replacement, char *newName);
void get_planarFrame(char *aa, long pnum, char *planar_atoms[], char *BDIR, char *idmsg,
                     long ib, long ie, char **AtomName, double **xyz, double *orien_i,
                     double *org_i);
void get_argFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);
void get_pheFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);
void get_tyrFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);
void get_trpFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);
void get_hisFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);
void get_lysFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);
void get_asnFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);
void get_glnFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);
void get_aspFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);
void get_gluFrame(char *BDIR, char *idmsg, long ib, long ie, char **AtomName,
                  double **xyz, double *orien_i, double *org_i);
void get_planarAA_frames(long num_residue, long **seidx, char **AtomName, char **ResName,
                         char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                         long *res_type, long *cidx);
void planarAA_blks(long num_residue, long **seidx, char **ResName, char *BDIR, long *cidx,
                   double **paa_orien, double **paa_org, char *alcfile);
void print_resid(long num_residue, long **seidx, char **ResName, char *ChainID,
                 long *ResSeq, char **Miscs, long *res_type);
void snap_atype(char **AtomName, long num_residue, long **seidx, long *res_type, long **atom_cidx);
long number_of_aa(long num_residue, long *res_type);
long number_of_nt(long num_residue, long *res_type);
void get_snap_par(double **rot1, double *org1, double **rot2, double *org2,
                  char *direction, double *trs_dist, double *rot_dist, double *pars,
                  double *orgP, double **rotP);
void write_snap_par(FILE * fp, char *direction, double dist, double rot_ang,
                    double *pars, double *orgP, double **rotP);
void write_atom_xyz(FILE * fp, char *fmt, double *xyz, double dft);
long set2frame(long inum, long *ivec, long **seidx, double **xyz, double *morg,
               double **mst, long *serial, double **xyz_pair);
long set2frameCa(long inum, long *ivec, long **seidx, long *res_type,
                 char **AtomName, double **xyz, double *morg, double **mst,
                 long *serial, double **xyz_pair);
long get_nextModelNumber(char *pdbfile);
void output_naa_str(char *filename, char *fmode, char *idmsg, long num, long *serial,
                    char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                    double **xyz_pair, char **Miscs, long out_org);
void cleanup_files(long renew, long cleanup);

/* cehs.c */

/* cmn_fncs.c */
long set_3letter_base_pdb(char *res_name, char *spdb);
void set_std_base_pdb(char *bdir, long irna, char bname, char *spdb);
void set_std_base_pdb00(char *bdir, long irna, char bname, char *spdb);
void print_used_time(time_t time0);
void parcat(char *str, double par, char *format, char *bstr);
void print_bp_crit(miscPars * misc_pars, FILE * fp);
char *my_getline(FILE * fp);
long csplit(char *str, char *item[], long itemsize, char sepc);
char *trim(char *a);
char *ltrim(char *a);
char *rtrim(char *a);
long itemize(char *str, char *item[], long itemsize);
long item_list(char *str, char *item[], long itemsize, char *sep_chars);
void refs_right_left(long bnum, double **orien, double **org, double **r1, double *o1,
                     double **r2, double *o2);
void refs_i_j(long b1, long b2, double *bp_orien, double *bp_org, double **r1,
              double *o1, double **r2, double *o2);
void ref_frame_i(long bnum, double *bp_orien, double *bp_org, double **r, double *o);
void mst2orien(double *orien_vec, long ioffset, double **mst);
void orien2mst(double *orien_vec, long ioffset, double **mst);
void x_y_z_2_mtx(double *x, double *y, double *z, double **mtx);
void mtx_2_x_y_z(double **mtx, double *x, double *y, double *z);
void cehs_average(long inum_base, long *ivec, double **orien, double **org, double **mst,
                  double *morg);
void geom_average(long inum_base, long *ivec, double **orien, double **org, double **mst,
                  double *morg);
void pair2mst(long inum_base, long *ivec, char **AtomName, char **ResName, char *ChainID,
              long *ResSeq, char **Miscs, double **xyz, double **orien, double **org,
              long **seidx, double *mst_orien, double *mst_org, long **htm_water,
              miscPars * misc_pars, FILE * fp);
void get_chi_angle(long num_residue, long *RY, char *bseq, long **seidx, double **xyz,
                   char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                   char **Miscs, double *chi, long **idxCN);
FILE *open_tmpfile(void);
FILE *open_file(char *filename, char *filemode);
long close_file(FILE * fp);
long exist_file(char *filename);
void remove_file(char *filename);
void rename_file(char *src, char *dst);
void copy_file_pointer(FILE * fpi, FILE * fpo, char *msg);
void cpcat_file(char *src, char *dst, char *method);
long upperstr(char *a);
long lowerstr(char *a);
char *my_strdup(const char *src);
void print_sep(FILE * fp, char x, long n);
void check_slash(char *BDIR);
void delete_end_slash(char *str);
char *basename(char *str);
long lround(double d);
void del_extension(char *fullname, char *okname);
void bname_noext(char *src, char *dst);
void fatal(char *fmt, ...);
void print_pdb_title(char *pdbfile, char *chain_list, FILE * fp);
long number_of_atoms(char *pdbfile, long hetatm, char *ALT_LIST);
long read_pdb(char *pdbfile, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs, long hetatm, char *ALT_LIST);
void free_pdb(long num, long *AtomSNum, char **AtomName, char **ResName, char *ChainID,
              long *ResSeq, double **xyz, char **Miscs);
void reset_xyz(long num, double **xyz, char *fmt);
void deduce_misc(char **Miscs, char **AtomName, long i, char *str);
long is_dna_with_backbone(long ib, long ie, char **AtomName);
void normalize_resName_atomName(long is_dna, const char *rname0, const char *aname0,
                                char *rname, char *aname);
void pdb_record(long ib, long ie, long *inum, long idx, char **AtomName, char **ResName,
                char *ChainID, long *ResSeq, double **xyz, char **Miscs, FILE * fp);
void write_pdb(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               double **xyz, char **Miscs, char *pdbfile);
void write_pdbml(long xml, long num, char **AtomName, char **ResName, char *ChainID,
                 long *ResSeq, double **xyz, char *pdbfile);
void write_pdbcnt(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                  double **xyz, long **connect, char *pdbfile);
void move_position(double **d, long nr, long nc, double *mpos);
long **residue_idx(long num, long *ResSeq, char **Miscs, char *ChainID, char **ResName,
                   long *num_residue);
long frag_contain_metal(long ib, long ie, long *is_metal);
void atom_metal(long num_atoms, char **AtomName, long *is_metal);
void residue_wtype(long num_residue, long **seidx, char **ResName, char **AtomName,
                   double **xyz, char **Miscs, long *res_type, long only_ntaa);
long residue_ident(char **AtomName, double **xyz, char **Miscs, long ib, long ie);
void normalize_atom_symbol(char *asym);
void get_atomlist(char **atomlist, long *num_sa);
long has_atom_name(long ib, long ie, char **AtomName, char *aname);
void get_baselist(char **baselist, long *num_sb);
void get_seq(long num_residue, long **seidx, char **AtomName, char **ResName,
             char *ChainID, long *ResSeq, char **Miscs, double **xyz, char *bseq, long *RY);
void get_bpseq(long ds, long num_bp, long **pair_num, long **seidx, char **AtomName,
               char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
               char **bp_seq, long *RY);
long strmatch_idx(char *str, char **strmat, long nb, long ne);
long num_strmatch(char *str, char **strmat, long nb, long ne);
void get_idmsg(char *rname, char cid, long snum, char icode, char *idmsg);
long find_1st_atom(char *str, char **strmat, long nb, long ne, char *idmsg);
double torsion(double **d);
double torsion2(double **d);
void get_BDIR(char *BDIR, char *filename);
void align2zaxis(long num, double *haxis, double **rotmat, double **xyz, double **xyzH);
void cov_matrix(double **a, double **b, long nr, long nc, double **cmtx);
double ls_fitting(double **sxyz, double **exyz, long n, double **fitted_xyz, double **R,
                  double *orgi);
void ls_plane(double **bxyz, long n, double *pnormal, double *ppos, double *odist, double *adist);
void arb_rotation(double *va, double ang_deg, double **rot_mtx);
double vec_ang(double *va, double *vb, double *vref);
void get_vector(double *va, double *vref, double deg_ang, double *vo);
void rotate(double **a, long i, long j, long k, long l, double *g, double *h, double s, double tau);
void eigsrt(double *d, double **v, long n);
void jacobi(double **a, long n, double *d, double **v);
void dludcmp(double **a, long n, long *indx, double *d);
void dlubksb(double **a, long n, long *indx, double *b);
void dinverse(double **a, long n, double **y);
void rotx(double ang_deg, double **rotmat);
void roty(double ang_deg, double **rotmat);
void rotz(double ang_deg, double **rotmat);
void get_alc_nums(char *alcname, long *num, long *nbond);
void read_alc(char *alcname, long *num, long *nbond, char **AtomName, double **xyz,
              long *ibase, long **linkage);
void write_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
               long **linkage, char *alcfile);
void free_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
              long zero_1, long **linkage);
void cnct_org(long num_bp, long ia, long ib, char **tAtomName, double **txyz,
              long *tibase, long **tlinkage, double **org_xyz);
void dsort(long n, double *a, long *idx);
void lsort(long n, long *a, long *idx);
void lreverse(long ia, long n, long *lvec);
void fig_title(FILE * fp);
void ps_title_cmds(FILE * fp, char *imgfile, long *bbox);
void get_fig_xy(long num, double **xyz, long nO, double **oxyz, long *urxy,
                long frame_box, FILE * fp);
void get_pxy(double *xy1, double *xy2, double r, double *px, double *py);
void alc2fig(long nobj, long *idx, long *depth, long **allobj, double **blkxyz,
             double **oxyz, long *ibase, long faces[][5], long *opts, FILE * fp);
void get_ps_xy(char *imgfile, long *urxy, long frame_box, FILE * fp);
void alc2ps(long nobj, long *idx, long **allobj, double **blkxyz, double **oxyz,
            long *ibase, long faces[][5], long *opts, FILE * fp);
void bring_atoms(long ib, long ie, long ra_num, char **AtomName, long *nmatch, long *batom);
void all_bring_atoms(long num_residue, long *RY, long **seidx, char **AtomName,
                     long *num_ring, long **ring_atom);
void base_idx(long num, char *bseq, long *ibase, long single);
long basepair_idx(char *bpi);
void plane_xyz(long num, double **xyz, double *ppos, double *nml, double **nxyz);
void prj2plane(long num, long ra_num, char **AtomName, double **xyz, double z0, double **nxyz);
void adjust_xy(long num, double **xyz, long nO, double **oxyz, double scale_factor,
               long default_size, long *urxy);
void get_depth(long nobj, long *zval, long *depth);
void raster3d_header(long num, double **xyz, double scale_factor, long no_header,
                     long frame_box, FILE * fp);
void get_r3dpars(double **base_col, double *hb_col, double *width3, double **atom_col,
                 char *label_style);
void r3d_rod(long itype, double *xyz1, double *xyz2, double rad, double *rgbv, FILE * fp);
void r3d_dash(double *xyz1, double *xyz2, double hb_width, double *hb_col, FILE * fp);
void r3d_sphere(double *xyz1, double rad, double *rgbv, FILE * fp);
void cpk_model(long num, long *idx, double **xyz, double ballrad, double **colrgb, FILE * fp);
void r3d_tripln(long itype, double *xyz1, double *xyz2, double *xyz3, double *rgbv, FILE * fp);
void r3d_block_edge(double *rgbv, long ioffset8, double **blkxyz, double w1, FILE * fp);
void base_label(double **rxyz, char *label_style, double *rgbv, char *bname_num, FILE * fp);
void fill_base_ring(long num_residue, long num_ring, long **ring_atom, double **xyz,
                    long *ibase, char *bseq, double **base_col, char *label_style,
                    long label_ring, long *ResSeq, FILE * fp);
void process_alc(char *alcfile, char *imgfile, double scale_factor, long *opts);
void alc_3images(long *opts, long nobj, long num_blk, long num_blk8, long nO_lkg,
                 long nO, double **oxyz, double **blkxyz, long *blkibase, long **linkage,
                 double scale_factor, char *imgfile);
void get_alc_objs(long num_blk, double **blkxyz, long nO, double **oxyz, long nO_lkg,
                  long **linkage, long faces[][5], long **allobj);
void get_fig_pars(double *dot_sep, long *dlcol, long *dwidth, long *bp1width,
                  long *bp2width, long **bc_idx, double *msat, double *Msat,
                  long *o_sides, long *line_width, long *join_style, long *cap_style, long *mfcol);
void frame_xyz(long side_view, double *morg, double **mst, long num, double **xyz);
void change_xyz(long side_view, double *morg, double **mst, long num, double **xyz);
void get_side_view(long ib, long ie, double **xyz);
void get_CNidx(long ds, long num_bp, long **chi, long **idx, long *nvec, long *C1b, long *C1e);
void add_3axes(long *num, char **AtomName, long *ibase, double **xyz, long *nbond,
               long **linkage, long side_view, double axis_len);
char *get_sequence(char *Wbase, long *num_bp);
char **single2double(long nbp, char *bseq, char *Wbase, char *Cbase);
long is_valid_base(char c, char *valid_bases);
long repeat_num(void);
char *read_sequence(char *seqfile, char *valid_bases, long *nbp);
char *read_repeat(char *crepeat, long fixed, char *valid_bases, long *nbp);
void combine_pstnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                    char **tAtomName, char **tResName, char *tChainID, long *tResSeq,
                    double **txyz, char **tAtomName2, double **txyz2);
void reverse_stnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                   char **tAtomName, char **tResName, char *tChainID, long *tResSeq,
                   double **txyz, char **tAtomName2, double **txyz2, long basep);
void pair_checking(long ip, long ds, long num_residue, char *pdbfile, long *num_bp,
                   long **pair_num);
void double_print_msg(char *msg, FILE * fp);
void drct_checking(long ds, long num_bp, long **pair_num, long **seidx, char **AtomName,
                   double **xyz, long *parallel, long *bbexist, long **o3p_brk, FILE * fp);
void residue_idstr(char chain_id, long res_seq, char *rname, char *idmsg);
void base_str(char chain_id, long res_seq, char *misc, char *rname, char bcode,
              long stnd, char *idmsg);
void write_lkglist(long nbond, long **linkage, char **AtomName, char **ResName,
                   char *ChainID, long *ResSeq, char **Miscs);
void hbond_info(long **pair_num, char *bseq, long **seidx, long *idx, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                long *RY, long *num_hbond, long **hb_linkage);
void hbond_pdb(long num, long num_residue, char *bseq, long **seidx, long *idx, long *RY,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, long *num_hbond, long **hb_linkage, long pwise);
void hbond_list(long i, long j, char **AtomName, char **ResName, char *ChainID,
                long *ResSeq, double **xyz, char **Miscs, char *bseq, long **seidx,
                long *idx, long **hb_linkage, miscPars * misc_pars, long **num_list,
                long *num_hbond, long ilayer, FILE * fph);
void hb_numlist(long i, long j, char basei, char basej, long **seidx, long *idx,
                char **AtomName, double **xyz, miscPars * misc_pars, long *num_hb, long **num_list);
void hb_information(long num_bp, long **pair_num, char **bp_seq, long **seidx, long *idx,
                    char **AtomName, double **xyz, long *WC_info, FILE * fp);
long good_hbatoms(miscPars * misc_pars, char *atom1, char *atom2, long idx1, long idx2);
void read_lkginfo(char *lkgfile, long num, long *nbond, long **linkage);
void read_hbinfo(char *hbfile, long num, long *num_hbond, long **hb_linkage);
void update_hb_idx(long idx, double *dtmp, long *ddidx, double *hb_dist, long cur_idx);
void hb_atompair(long num_hbonds, char **hb_atom1, char **hb_atom2, double *hb_dist,
                 long *lkg_type, miscPars * misc_pars);
long validate_hbonds(long num_hbonds, double *hb_dist, long *lkg_type, char *hb_type,
                     char basei, char basej, char **hb_atom1, char **hb_atom2);
void get_hbond_ij(long i, long j, char basei, char basej, miscPars * misc_pars,
                  long **seidx, long *idx, char **AtomName, double **xyz, char *hb_info);
char donor_acceptor(char basei, char basej, char *hb_atom1, char *hb_atom2);
long asym_idx(char *asym, char atoms_list[NELE][3], long dft_lval);
void atom_info(long idx, char atoms_list[NELE][3], double *covalence_radii, double *vdw_radii);
void aname2asym(const char *aname0, char *my_asym, long num_sa, char **atomlist);
void atom_idx(long num, char **AtomName, char **Miscs, long *idx);
void get_bonds(long num, char **AtomName, double **xyz, long num_residue, long *RY,
               long **seidx, long **connect);
void atom_linkage(long ib, long ie, long *idx, double **xyz, char **Miscs, char *ChainID,
                  long nbond_estimated, long *nbond, long **linkage);
void lkg2connect(char **AtomName, long ib, long ie, long nbond, long **linkage, long **connect);
void init_htm_water(long waters, long num, long num_residue, long *idx, long **htm_water);
void identify_htw(long num_residue, long **seidx, long *RY, char **AtomName,
                  char **ResName, char *ChainID, long *ResSeq, char **Miscs,
                  double **xyz, long **htm_water);
long attached_residues(long inum_base, long *ivec, long *ivec2, long **seidx,
                       double **xyz, long **htm_water, miscPars * misc_pars);
void print_pairinfo(long i, long j, char basei, char basej, double *rtn_val, double *chi,
                    miscPars * misc_pars, long **seidx, long *idx, char **AtomName,
                    double **xyz, char *bseq, long detailed, FILE * fp);
void check_pair(long i, long j, char *bseq, long **seidx, double **xyz,
                double **NC1xyz, double **orien, double **org, long *idx,
                char **AtomName, miscPars * misc_pars, double *rtn_val,
                long *bpid, long **ring_atom, long network);
void o3_p_xyz(long ib, long ie, char *aname, char **AtomName, double **xyz,
              double *o3_or_p, long idx);
long is_baseatom(char *atomname);
void base_info(long num_residue, char *bseq, long **seidx, long *RY,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, double **orien, double **org,
               double **NC1xyz, double **o3_p);
void help3dna_usage(char *program_name);
void help3dna(char *program_name);
void delH_pdbfile(char *inpfile, char *outfile);
void contact_msg(long prt_msg);
long str_pmatch(char *str, char *sstr);
long case_str_pmatch(char *str, char *sstr);
long is_numeric(char *str);
double cvt2double(char *str);
long cvt2long(char *str);
long equalsign_pos(char *str);
long get_lvalue(char *str, long vmin, long vmax);
double get_dvalue(char *str, double vmin, double vmax);
void get_strvalue(char *str, char *dst, long expand_tilde);
void reverse_string(char *str);
void cvtstr_set1toc2(char *str, char *set1, char c2);
void cvtstr_c1toc2(char *str, char c1, char c2);
double z1_z2_angle_in_0_to_90(double *z1, double *z2);
void do_nothing(void);
void skip_lines(long num, FILE * fp);
void check_havefile(char *filename, char *msg);

/* comb_str.c */

/* ex_str.c */

/* fiber.c */

/* find_pair.c */

/* fncs_slre.c */
const char *slre_match(enum slre_option options, const char *re, const char *buf, int buf_len, ...);
int lux_match(enum slre_option options, const char *re, const char *buf);
int lux_bcmatch(const char *buf, const char *re);
int lux_ncmatch(const char *buf, const char *re);

/* frame_mol.c */

/* get_part.c */

/* mutate_bases.c */

/* nrutil.c */
void nrerror(char *error_text);
void vector_boundary_check(long nl, long nh, char *fun_name);
void matrix_boundary_check(long nrl, long nrh, long ncl, long nch, char *fun_name);
char *cvector(long nl, long nh);
char *cvector_nr(long nl, long nh);
double *dvector(long nl, long nh);
double *dvector_nr(long nl, long nh);
long *lvector(long nl, long nh);
long *lvector_nr(long nl, long nh);
char **cmatrix(long nrl, long nrh, long ncl, long nch);
char **cmatrix_nr(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix_nr(long nrl, long nrh, long ncl, long nch);
long **lmatrix(long nrl, long nrh, long ncl, long nch);
long **lmatrix_nr(long nrl, long nrh, long ncl, long nch);
void free_cvector(char *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_lvector(long *v, long nl, long nh);
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch);
double dval_sqr(double dval);
void dval_swap(double *pa, double *pb);
void lval_swap(long *pa, long *pb);
void cval_swap(char *pa, char *pb);
double dval_max(double a, double b);
double dval_min(double a, double b);
long lval_max(long a, long b);
long lval_min(long a, long b);
double abs_dval_diff(double a, double b);
long lval_in_set(long lval, long ib, long ie, long *s);
long dval_in_range(double dval, double dlow, double dhigh);
long lval_in_range(long lval, long llow, long lhigh);
void max_dmatrix(double **d, long nr, long nc, double *maxdm);
void min_dmatrix(double **d, long nr, long nc, double *mindm);
void ave_dmatrix(double **d, long nr, long nc, double *avedm);
void std_dmatrix(double **d, long nr, long nc, double *stddm);
double max_dvector(double *d, long nl, long nh);
double min_dvector(double *d, long nl, long nh);
double ave_dvector(double *d, long n);
double std_dvector(double *d, long n);
void init_cmatrix(char **cmtx, long nrl, long nrh, long ncl, long nch, char init_val);
void init_dmatrix(double **dmtx, long nrl, long nrh, long ncl, long nch, double init_val);
void init_lmatrix(long **lmtx, long nrl, long nrh, long ncl, long nch, long init_val);
void init_cvector(char *cvec, long ib, long ie, char init_val);
void init_cvector_all(char *cvec, long ib, long ie, char init_val);
void init_dvector(double *dvec, long ib, long ie, double init_val);
void init_lvector(long *lvec, long ib, long ie, long init_val);
void copy_dvector(double *d, double *s, long nl, long nh);
void copy_lvector(long *d, long *s, long nl, long nh);
int dval_compare(const void *v1, const void *v2);
int lval_compare(const void *v1, const void *v2);
int cstr_compare(const void *v1, const void *v2);
void negate_xyz(double *xyz1);
double p1p2_dist(double *xyz1, double *xyz2);
void p1p2_ave(double *xyz1, double *xyz2, double *ave);
long within_limits(double *xyz1, double *xyz2, double dlow, double dhigh);
void sumxyz(double *xyz1, double *xyz2, double *sxyz);
void avexyz(double *xyz1, double *xyz2, double *mxyz);
void ddxyz(double *xyz1, double *xyz2, double *dxyz);
void cpxyz(double *xyz1, double *xyz2);
void vec_orth(double *va, double *vref);
double dot(double *va, double *vb);
void cross(double *va, double *vb, double *vc);
long sign_control(double *va, double *vb, double *vref);
double veclen(double *va);
void vec_norm(double *va);
double dot2ang(double dotval);
double magang(double *va, double *vb);
double rad2deg(double ang);
double deg2rad(double ang);
void copy_dmatrix(double **a, long nr, long nc, double **o);
void copy_lmatrix(long **a, long nr, long nc, long **o);
void multi_matrix(double **a, long nra, long nca, double **b, long nrb, long ncb, double **o);
void multi_vec_matrix(double *a, long n, double **b, long nr, long nc, double *o);
void multi_vec_Tmatrix(double *a, long n, double **b, long nr, long nc, double *o);
void transpose_matrix(double **a, long nr, long nc, double **o);
void identity_matrix(double **d, long n);

/* o1p_o2p.c */

/* pdb2img.c */

/* r3d_atom.c */

/* reb_fncs.c */
void link_o3_p(long num_residue, long **seidx, char **AtomName, double **xyz,
               char *ChainID, long **connect);
void atom_lkg(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
              double **xyz, long xml, char *outfile);
void atomic_pdb1(long num_bp, long num_atoms, long num_max_per_residue, long is_helical,
                 long xdir, char **bp_seq, double **step_par, char *BDIR, long xml, char *outfile);
void atomic_pdb2(long parallel, long num_bp, long num_atoms, long num_max_per_residue,
                 long is_helical, long xdir, char **bp_seq, double **bp_par,
                 double **step_par, char *BDIR, long xml, char *outfile);
void base_c1_atoms(char **AtomName, long ib, long ie, long *num_batom, long *batom);
void extract_base_atoms(long num, char **AtomName, double **xyz, long *bnum,
                        char **bAtomName, double **bxyz);
void atomic_base_p(long parallel, long num_bp, long num_atoms, long num_max_per_residue,
                   long is_helical, long xdir, char **bp_seq, double **bp_par,
                   double **step_par, char *BDIR, long *pidx, long xml, char *outfile);
void num_PDB_atoms(long num_bp, long is_single, char **bp_seq, char *BDIR,
                   long *num_atoms, long *num_max_per_residue);
void set_bp_pdb(long num_max_per_residue, char *fname1, char *fname2, double *param,
                long xdir, long *num1, long *num2, char **AtomName1, char **AtomName2,
                double **xyz1, double **xyz2, char ap);
void base_fname(char bname, char *BDIR, char *fname);
void block_alc1(long num_bp, long is_single, long is_helical, long xdir, char **bp_seq,
                double **step_par, char *BDIR, char *outfile);
void block_alc2(long num_bp, long is_helical, long xdir, char **bp_seq, double **bp_par,
                double **step_par, char *BDIR, char *outfile);
void set_bp_alc(char *fname1, char *fname2, double *param, long xdir, long num,
                long nbond, double **bp_xyz, char ap);
void xbpfunc(double *param, double **orien, double **mst, double *pos, double *mpos);
void xhelfunc(double *param, double **orien, double **mst, double *pos, double *mpos);
void print_ref_frames(FILE * fp, double *pos_next, double **orien_next);

/* rebuild.c */

/* regular_dna.c */

/* rotate_mol.c */

/* stack2img.c */

/* std_base.c */

/* step_hel.c */
