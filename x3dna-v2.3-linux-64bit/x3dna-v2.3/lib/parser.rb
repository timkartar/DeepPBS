### 3DNA parser for analysis results of a MODEL/ENDMDL delineated ensemble
###     3DNA v2.3-2017feb08, created and maintained by Xiang-Jun Lu (PhD)

$VERBOSE = true
abort "Please set the 3DNA environment variable" unless ENV['X3DNA']
require ENV['X3DNA'] + '/lib/miscs'

class Parse_3DNA
    def extract_pars(parlist, parmtx, pars)
        parlist.each_with_index do |p, i|
            pars[p] = parmtx.map { |x| x[i] }
        end
    end

    def parse_hbond_parameters(aFile, num_bp, pars)
        parmtx = []
        num_bp.times { parmtx << aFile.gets.chomp[14..-1] } # [2]  N6 - O4  3.01  N1 - N3  2.83
        pars['hbond'] = parmtx
    end

    def parse_origin_normalVector(aFile, num_bp, pars, par_names)
        parmtx = []
        Utils.skip_lines(aFile, 3)
        num_bp.times do
            str = aFile.gets.split
            parmtx << [str[2..4].join(" "), str[5..7].join(" ")]
        end
        extract_pars(par_names, parmtx, pars)
    end

    def parse_overlap_area(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 16)
        (num_bp - 1).times do
            str = aFile.gets
            parmtx << [str[63, 5], str[69, 5]]
        end
        extract_pars(%w(area area0), parmtx, pars)
    end

    def parse_overlap_area_single(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 16)
        (num_bp - 1).times do
            str = aFile.gets
            parmtx << [str[11, 5], str[17, 5]]
        end
        extract_pars(%w(area area0), parmtx, pars)
    end

    def parse_base_pair_parameters(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile)
        num_bp.times { parmtx << aFile.gets.split[2, 6] }
        extract_pars(%w(shear stretch stagger buckle propeller opening), parmtx, pars)
    end

    def parse_step_parameters(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile)
        (num_bp - 1).times { parmtx << aFile.gets.split[2, 6] }
        extract_pars(%w(shift slide rise tilt roll twist), parmtx, pars)
    end

    def parse_helical_parameters(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile)
        (num_bp - 1).times { parmtx << aFile.gets.split[2, 6] }
        extract_pars(%w(x_displacement y_displacement h_rise inclination tip h_twist), parmtx, pars)
    end

    def parse_lambda_NC_dists(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 7)
        num_bp.times { parmtx << aFile.gets.split[2, 5] }
        extract_pars(%w(lambda1 lambda2 dist_c1_c1 dist_rn9_yn1 dist_rc8_yc6), parmtx, pars)
    end

    def parse_zp_zph(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 3)
        (num_bp - 1).times { parmtx << aFile.gets.split.values_at(4, 7) }
        extract_pars(%w(zp zph), parmtx, pars)
    end

    def parse_groove_widths(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 10)
        (num_bp - 1).times { parmtx << aFile.gets.split[2, 4] }
        extract_pars(%w(minor_gw_pp minor_gw_refined major_gw_pp major_gw_refined), parmtx, pars)
    end

    def parse_torsion_angles(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 13)
        num_bp.times { parmtx << aFile.gets.split[2, 7] }
        extract_pars(%w(alpha1 beta1 gamma1 delta1 epsilon1 zeta1 chi1), parmtx, pars)

        parmtx = []     # for strand II, 3'--->5' direction
        Utils.skip_lines(aFile, 3)
        num_bp.times { parmtx << aFile.gets.split[2, 7] }
        extract_pars(%w(alpha2 beta2 gamma2 delta2 epsilon2 zeta2 chi2), parmtx, pars)
    end

    def parse_torsion_angles_single(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 12)
        num_bp.times { parmtx << aFile.gets.split[2, 7] }
        extract_pars(%w(alpha beta gamma delta epsilon zeta chi), parmtx, pars)
    end

    def parse_sugar_puckers(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 12)
        num_bp.times { parmtx << aFile.gets.split[7, 3] }
        extract_pars(%w(tm1 phase1 puckering1), parmtx, pars)

        parmtx = []     # for strand II, 3'--->5' direction
        Utils.skip_lines(aFile, 3)
        num_bp.times { parmtx << aFile.gets.split[7, 3] }
        extract_pars(%w(tm2 phase2 puckering2), parmtx, pars)
    end

    def parse_sugar_puckers_single(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 11)
        num_bp.times { parmtx << aFile.gets.split[7, 3] }
        extract_pars(%w(tm phase puckering), parmtx, pars)
    end

    def parse_same_chain_PC_dists(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 3)
        (num_bp - 1).times { parmtx << aFile.gets.split.values_at(2, 3, 6, 7) }
        extract_pars(%w(dist_pp1 dist_cc1 dist_pp2 dist_cc2), parmtx, pars)
    end

    def parse_same_chain_PC_dists_single(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 3)
        (num_bp - 1).times { parmtx << aFile.gets.split.values_at(2, 3) }
        extract_pars(%w(dist_pp dist_cc), parmtx, pars)
    end

    def parse_helix_radius(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 4)
        (num_bp - 1).times { parmtx << aFile.gets.split[2..7] }
        extract_pars(%w(radius1_p radius1_o4 radius1_c1
                        radius2_p radius2_o4 radius2_c1), parmtx, pars)
    end

    def parse_helix_radius_single(aFile, num_bp, pars)
        parmtx = []
        Utils.skip_lines(aFile, 3)
        (num_bp - 1).times { parmtx << aFile.gets.split[2..4] }
        extract_pars(%w(radius_p radius_o4 radius_c1), parmtx, pars)
    end

    def parse_helix_pos_vector(aFile, num_bp, pars, num=3)
        parmtx = []
        Utils.skip_lines(aFile, num)
        (num_bp - 1).times do
            str = aFile.gets.split
            parmtx << [str[2..4].join(" "), str[5..7].join(" ")]
        end
        extract_pars(%w(hx_pos hx_vec), parmtx, pars)
    end

    def parse_3dna_output
        pars = {} # key: parameter (e.g. 'roll'); value: array of steps/nucleotides

        File.open(@parfile) do |aFile|
            while line = aFile.gets do
                if line =~ /Detailed H-bond information/
                    parse_hbond_parameters(aFile, @num_bp, pars)
                elsif line =~ /Overlap area in Angstrom/
                    parse_overlap_area(aFile, @num_bp, pars)
                elsif line =~ /Origin /
                    parse_origin_normalVector(aFile, @num_bp, pars, %w(bpo_xyz bpn_vec))
                elsif line =~ /Local base-pair parameters/
                    parse_base_pair_parameters(aFile, @num_bp, pars)
                elsif line =~ /Local base-pair step parameters/
                    parse_step_parameters(aFile, @num_bp, pars)
                elsif line =~ /Local base-pair helical parameters/
                    parse_helical_parameters(aFile, @num_bp, pars)
                elsif line =~ /Strand I base ring /
                    parse_origin_normalVector(aFile, @num_bp, pars, %w(b1_cxyz b1_nvec))
                elsif line =~ /Strand II base ring /
                    parse_origin_normalVector(aFile, @num_bp, pars, %w(b2_cxyz b2_nvec))
                elsif line =~ /lambda: virtual angle between/
                    parse_lambda_NC_dists(aFile, @num_bp, pars)
                elsif line =~ /Classification of each dinucleotide step/
                    parse_zp_zph(aFile, @num_bp, pars)
                elsif line =~ /Minor and major groove widths/
                    parse_groove_widths(aFile, @num_bp, pars)
                elsif line =~ /Main chain and chi torsion angles/
                    parse_torsion_angles(aFile, @num_bp, pars)
                elsif line =~ /Sugar conformational parameters/
                    parse_sugar_puckers(aFile, @num_bp, pars)
                elsif line =~ /Same strand P--P/
                    parse_same_chain_PC_dists(aFile, @num_bp, pars)
                elsif line =~ /Helix radius /
                    parse_helix_radius(aFile, @num_bp, pars)
                elsif line =~ /Position /
                    parse_helix_pos_vector(aFile, @num_bp, pars)
                    break           # no more ...
                end
            end
        end

        pars
    end

    def initialize(num_bp, parfile)
        @num_bp = num_bp

        Utils.fatal("File #{parfile} nonexistent") unless File.exist?(parfile)
        Utils.fatal("File #{parfile} empty") if File.zero?(parfile)
        @parfile = parfile
    end

    ### ******************* for option -single *******************
    def parse_3dna_single_output
        pars = {} # key: parameter (e.g. 'chi'); value: array of steps/nucleotides
        num_nt = @num_bp

        File.open(@parfile) do |aFile|
            while line = aFile.gets do
                if line =~ /Overlap area in Angstrom/
                    parse_overlap_area_single(aFile, num_nt, pars)
                elsif line =~ /Local base step parameters/
                    parse_step_parameters(aFile, num_nt, pars)
                elsif line =~ /Local base helical parameters/
                    parse_helical_parameters(aFile, num_nt, pars)
                elsif line =~ /Base ring center /
                    parse_origin_normalVector(aFile, num_nt, pars, %w(b_cxyz b_nvec))
                elsif line =~ /Main chain and chi torsion angles/
                    parse_torsion_angles_single(aFile, num_nt, pars)
                elsif line =~ /Sugar conformational parameters/
                    parse_sugar_puckers_single(aFile, num_nt, pars)
                elsif line =~ /Same strand P--P/
                    parse_same_chain_PC_dists_single(aFile, num_nt, pars)
                elsif line =~ /Helix radius /
                    parse_helix_radius_single(aFile, num_nt, pars)
                elsif line =~ /Position /
                    parse_helix_pos_vector(aFile, num_nt, pars, 2)
                    break           # no more ...
                end
            end
        end

        pars
    end

    ### ******************* for option -torsion *******************
    def parse_torsion_BI_BII(aFile, num_nt, pars)
        parmtx = []
        Utils.skip_lines(aFile, 19)
        num_nt.times { parmtx << aFile.gets.split[2, 10] }
        extract_pars(%w(chi chi_type alpha beta gamma delta epsilon zeta epsilon_zeta bb_type),
                     parmtx, pars)
    end

    def parse_pseudo_torsion(aFile, num_nt, pars)
        parmtx = []
        Utils.skip_lines(aFile, 11)
        num_nt.times { parmtx << aFile.gets.split[2, 6] }
        extract_pars(%w(eta theta eta_c1 theta_c1 eta_baseorg theta_baseorg), parmtx, pars)
    end

    def parse_sugar_torsion_Zp(aFile, num_nt, pars)
        parmtx = []
        Utils.skip_lines(aFile, 18)
        num_nt.times { parmtx << aFile.gets.split[2, 10] }
        extract_pars(%w(v0 v1 v2 v3 v4 tm phase puckering zp_base dp_ncbond), parmtx, pars)
    end

    def parse_3dna_torsion_output
        pars = {} # key: parameter (e.g. 'chi'); value: array of nucleotides
        num_nt = @num_bp

        File.open(@parfile) do |aFile|
            while line = aFile.gets do
                if line =~ /Main chain and chi torsion angles/
                    parse_torsion_BI_BII(aFile, num_nt, pars)
                elsif line =~ /Pseudo /
                    parse_pseudo_torsion(aFile, num_nt, pars)
                elsif line =~ /Sugar conformational parameters/
                    parse_sugar_torsion_Zp(aFile, num_nt, pars)
                    break           # no more ...
                end
            end
        end

        pars
    end
end
