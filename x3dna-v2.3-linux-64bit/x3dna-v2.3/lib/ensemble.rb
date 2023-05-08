### Module and classes to handle MODEL/ENDMDL delineated ensembles
###     3DNA v2.3-2017feb08, created and maintained by Xiang-Jun Lu (PhD)

$VERBOSE = true
require 'fileutils'

abort "Please set the 3DNA environment variable" unless ENV['X3DNA']
require ENV['X3DNA'] + '/lib/trollop'
require ENV['X3DNA'] + '/lib/miscs'
require ENV['X3DNA'] + '/lib/parser'

module Ensemble
    def self.set_globals        # analyze/reorient
        $x3dna = Utils.x3dna_bindir

        $temp_model = "temp_model" # temporary model file (.pdb/.inp/.out)
        $model_list = "model_list.dat" # default models list
        $pdb_list = "pdb_list.dat"     # default PDB list
        $msg = "2> msgfile" # stderr redirected to 'msgfile'
        $oname = "ensemble_example" # basename of the output parameters file

        $errlog = false         # switch to log error or abort program
        $errcmd = false         # error in running command

        Dir["#{$temp_model}.*"].each { |t| File.delete(t) } # Utils.dcmnfile
    end

    class Analyze
        def main
            Ensemble.set_globals
            opts = parse_options

            $errlog = opts[:errlog]

            if opts[:ensemble]
                ensemble = opts[:ensemble]
                models = Utils.extract_model_numbers_from_ensemble_pdb(ensemble)
                Utils.fatal("Ensemble #{ensemble} contains #{models.size} models: " +
                            "#{models.first}...#{models.last}") if opts[:info]
            elsif opts[:models]
                ensemble, models = Utils.readlist_model_numbers(opts[:models])
            elsif opts[:list]
                pdblist = Utils.collect_pdb_list(opts[:list])
            elsif opts[:pattern] # "pdbdir/*.pdb" (quoted)
                pdblist = Dir[opts[:pattern]]
            else                # opts[:one]
                ensemble = opts[:one]
                models = [1]    # only one structure, taken as model #1
            end

            if models           # opts[:ensemble] or opts[:models]
                Utils.fatal("no MODEL/ENDMDL pairs in ensemble file") if models.size == 0
                Utils.output_model_list(ensemble, models, $model_list)
                process_ensemble_models(ensemble, models, opts)
            else                # opts[:list] or opts[:pattern]
                Utils.fatal("no matched PDB files") if pdblist.size == 0
                Utils.output_pdb_list(pdblist, $pdb_list)
                process_pdb_list(pdblist, opts)
            end
        end

        def write_out_parameters(all_pars, outfile, entries)
            File.open(outfile, "w") do |aFile|
                aFile.puts <<TXT
# Detailed output generated with '#{$prg} analyze' of a MODEL/ENDMDL
#     delineated ensemble of DNA/RNA structures. To extract parameters, run:
#     '#{$prg} extract'
TXT
                all_pars[0].keys.sort.each do |pname|
                    aFile.puts "\n<#{pname}>\t\t# with #{all_pars[0][pname].size} data columns"
                    all_pars.each_with_index do |pdata, idx| # loop over each model
                        next unless pdata
                        c0 = File.basename(entries[idx].to_s.downcase, ".*")
                        str = pdata[pname].map { |x| x.sub(/---/, 'na') }.join("\t")
                        aFile.puts "#{c0}\t#{str}"
                    end
                    aFile.puts "</#{pname}>"
                end
            end
        end

        def run_analyze_parse_parameters(opts, renew_bpfile=false)
            $num_bp = Utils.update_bpfile(opts[:bpfile], $temp_model) if renew_bpfile
            chk_ring = opts[:ring] ? "-ring" : ""
            Utils.run_cmd("#{$x3dna}analyze -c #{chk_ring} #{$temp_model}.inp #{$msg}")
            $errcmd ? nil : Parse_3DNA.new($num_bp, "#{$temp_model}.out").parse_3dna_output
        end

        def run_single_strand_analyze_parse_parameters(opts, renew_ntlist=false)
            $num_nt = Utils.update_ntlist($temp_model) if renew_ntlist
            chk_ring = opts[:ring] ? "-ring" : ""
            Utils.run_cmd("#{$x3dna}analyze #{chk_ring} #{$temp_model}.nts #{$msg}")
            $errcmd ? nil : Parse_3DNA.new($num_nt, "#{$temp_model}.outs").parse_3dna_single_output
        end

        def run_torsion_analyze_parse_parameters(opts, renew_ntlist=false)
            $num_nt = Utils.update_ntlist($temp_model) if renew_ntlist
            Utils.run_cmd("#{$x3dna}analyze -tor=#{$temp_model}.tor #{$temp_model}.pdb #{$msg}")
            $errcmd ? nil : Parse_3DNA.new($num_nt, "#{$temp_model}.tor").parse_3dna_torsion_output
        end

        # all_pars[] is an array, each row contains an hash for a model. The
        # hash keys are the extracted parameters, e.g., 'slide'. The value of
        # each corresponding key is an array, containing parameters of all
        # pairs/steps in a model.
        def process_ensemble_models(ensemble, models, opts)
            all_pars = []
            tnum = models.size

            models.each_with_index do |model, num|
                puts "Process model ##{model}\t#{num + 1} / #{tnum}"
                cmd = "#{$x3dna}ex_str -nmr -#{model} #{ensemble} #{$temp_model}.pdb"
                Utils.run_cmd("#{cmd} #{$msg}")
                all_pars <<
                    if opts[:single]
                        run_single_strand_analyze_parse_parameters(opts, num == 0)
                    elsif opts[:torsion]
                        run_torsion_analyze_parse_parameters(opts, num == 0)
                    else
                        run_analyze_parse_parameters(opts, num == 0)
                    end
            end

            write_out_parameters(all_pars, opts[:outfile], models)
        end

        def process_pdb_list(pdblist, opts)
            all_pars = []
            tnum = pdblist.size

            pdblist.each_with_index do |pdb, num|
                puts "Process PDB file #{pdb}\t#{num + 1} / #{tnum}"
                FileUtils.cp(pdb, "#{$temp_model}.pdb")
                all_pars <<
                    if opts[:single]
                        run_single_strand_analyze_parse_parameters(opts, num == 0)
                    elsif opts[:torsion]
                        run_torsion_analyze_parse_parameters(opts, num == 0)
                    else
                        run_analyze_parse_parameters(opts, num == 0)
                    end
            end

            write_out_parameters(all_pars, opts[:outfile], pdblist)
        end

        def docinfo
            sub = 'analyze'
            <<TXT
#{$sline}
Analyze a MODEL/ENDMDL delineated ensemble of NMR structures or MD
trajectories. All models must correspond to different conformations
of the same molecule. For the analysis of duplexes (default), a template
base-pair input file, generated with 'find_pair' and manually edited
as necessary, must be provided.

Usage:
        #{$prg} #{sub} options
Examples:
        #{$prg} #{sub} -b bpfile.dat -e sample_md0.pdb
             # 21 models (0-20); output (default): '#{$oname}.out'
             # also generate '#{$model_list}', see example below
        #{$prg} #{sub} -b bpfile.dat -m #{$model_list} -o #{$oname}2.out
             # diff #{$oname}.out #{$oname}2.out

        #{$prg} #{sub} -b bpfile.dat -p 'pdbdir/model_*.pdb' -o #{$oname}3.out
             # note to quote the -p option; 20 models (1-20)
             # also generate '#{$pdb_list}', see example below
        #{$prg} #{sub} -b bpfile.dat -l #{$pdb_list} -o #{$oname}4.out
             # diff #{$oname}3.out #{$oname}4.out
             # note the order of the models: 1, 10..19, 2, 20, 3..9

        #{$prg} #{sub} -s -e sample_md0.pdb
             # perform a 'single'-stranded analysis
        #{$prg} #{sub} -t -e sample_md0.pdb
             # calculate all 'torsion' angles

        find_pair 355d.pdb 355d.bps
        #{$prg} #{sub} -b 355d.bps --one 355d.pdb
             # process the structure file 355d.pdb specified in 355d.bps
Options:
#{$sline}
TXT
        end

        def parse_options
            txt = docinfo
            opts = Trollop::options do
                banner txt      # cannot directly call docinfo
                opt :bpfile, "Name of file containing base-pairing info", :type => :string
                opt :outfile, "Output file", :default => "#{$oname}.out"

                opt :single, "Single-stranded DNA/RNA", :default => false
                opt :torsion, "Torsion angles", :default => false
                opt :ring, "Base ring center & normal vector", :default => false

                opt :ensemble, "Ensemble delineated with MODEL/ENDMDL pairs", :type => :string
                opt :models, "File containing an explicit list of model numbers", :type => :string
                opt :pattern, "Pattern of model files to process (e.g., *.pdb)", :type => :string
                opt :list, "File containing an explicit list of models", :type => :string
                opt :one, "One regular structure [special case]", :type => :string

                opt :info, "Show only model info in the ensemble [with -e]", :default => false

                opt :errlog, "Print error message instead of abort program", :default => false
            end
            k = Utils.count_set_options([:ensemble, :models, :pattern, :list, :one], opts)
            Utils.fatal("Specify only one of: ensemble|models|pattern|list|one") unless k == 1

            Trollop::die :ensemble, "must exist" if opts[:ensemble] && !File.exist?(opts[:ensemble])
            Trollop::die :models, "must exist" if opts[:models] && !File.exist?(opts[:models])
            Trollop::die :list, "must exist" if opts[:list] && !File.exist?(opts[:list])
            Trollop::die :one, "must exist" if opts[:one] && !File.exist?(opts[:one])

            unless (opts[:info] && opts[:ensemble]) || opts[:single] || opts[:torsion]
                Trollop::die :bpfile, "must be specified" unless opts[:bpfile]
                Trollop::die :bpfile, "cannot be empty" if opts[:bpfile].empty?
                Trollop::die :bpfile, "must exist" unless File.exist?(opts[:bpfile])
            end

            opts
        end
    end

    class Extract
        def main
            Ensemble.set_globals
            opts = parse_options

            parslst = list_of_parameters(opts[:fromfile])
            clean_parfiles_and_exit(parslst) if opts[:clean]
            list_pars_and_exit(parslst) if opts[:list]

            if opts[:par_name]
                select_one_parameter(parslst, opts)
            else
                select_all_parameters(parslst, opts)
            end
        end

        def clean_parfiles_and_exit(parslst)
            parslst.each do |par|
                parfile = "#{$oname}_#{par}.out"
                File.delete(parfile) if File.exist?(parfile)
            end
            Utils.fatal("")
        end

        def list_pars_and_exit(parslst)
            puts $sline
            parslst.each_with_index do |par, idx|
                k = idx + 1
                str = "[%2d] %-20s" % [k, par]
                print str
                puts if k % 3 == 0
            end
            Utils.fatal("\n#{$sline}")
        end

        def select_one_parameter(parslst, opts)
            par = get_parameter_name(opts[:par_name], parslst)
            partxt = extract_one_parameter(par, opts)
            if opts[:outfile] == "stdout"
                puts partxt
            else
                File.open(opts[:outfile], "w") do |aFile|
                    aFile.puts partxt
                end
            end
        end

        def select_all_parameters(parslst, opts)
            parslst.each do |par|
                partxt = extract_one_parameter(par, opts)
                File.open("#{$oname}_#{par}.out", "w") do |aFile|
                    aFile.puts partxt
                end
            end
        end

        def get_parameter_name(par_name, parslst)
            max = parslst.size
            idx = par_name.to_i

            if idx > max || idx < 0
                Utils.fatal("Parameter index number '#{idx}' out of range: 1..#{max}")
            elsif (1..max).include?(idx)
                okidx = idx - 1
            else                # name
                matches = []
                parslst.each_with_index do |p, k|
                    matches << k if p =~ /^#{par_name}/i  # case insensitive
                end
                num = matches.size
                if num == 1
                    okidx = matches.first
                elsif num > 1
                    puts "Input '-p #{par_name}' matches the following #{num} parameters:"
                    matches.each do |k|
                        puts "\t[%2d] %-20s" % [k + 1, parslst[k]]
                    end
                    Utils.fatal("Please be specific and try again.")
                else
                    Utils.fatal("Input '-p #{par_name}' matches nothing.\n" +
                                "Please try: '#{$prg} -l' to see a full list of parameters.")
                end
            end

            par = parslst[okidx] # full parameter name
            $stderr.puts "You've pick paramter: [%2d] %-20s" % [okidx + 1, par]

            par
        end

        def list_of_parameters(parsfile)
            parslst = []
            File.open(parsfile).each_line do |line|
                parslst << $1.downcase if line =~ /^\s*<\/(.*)>/
            end
            Utils.fatal("No parameters found in <#{parfile}>") if parslst.empty?

            parslst.sort
        end

        def extract_one_parameter(par, opts)
            partxt = []
            found = false
            ib, ie = opts[:end_bps]

            File.open(opts[:fromfile]).each_line do |line|
                if line =~ /^\s*<#{par}>/i
                    found = true
                    next
                end
                break if line =~ /^\s*<\/#{par}>/i
                if found
                    items = line.split("\t")
                    if opts[:no_1col]
                        selected = items[(ib+1)..(-ie-1)]
                    else        # with first column for row label
                        selected = items.values_at(0, (ib+1)..(-ie-1))
                    end
                    partxt << selected
                end
            end

            Utils.fatal("Parameter <#{par}> not found in '#{opts[:fromfile]}'") unless found

            fmtdat = []
            partxt.each { |str| fmtdat << str.join(opts[:separator]) }

            fmtdat
        end

        def docinfo
            sub = 'extract'
            <<TXT
#{$sline}
Extract 3DNA structural parameters of an ensemble of NMR structures or
MD trajectories, after running '#{$prg} analyze'. The extracted
parameters are intended to be exported into Excel, Matlab and R etc for
further data analysis/visualization.

Usage:
        #{$prg} #{sub} options
Examples:
        #{$prg} #{sub} -l
             # to see a list of all parameters
        #{$prg} #{sub} -p prop
             # for propeller, no need to specify full: -p pr suffices
             # -p 36 also fine (see above); use '#{$oname}.out'
        #{$prg} #{sub} -p slide -s , -f #{$oname}3.out
             # comma separated, from file '#{$oname}3.out'
        #{$prg} #{sub} -p roll -s ' ' -n -o roll.dat
             # space separated, no row-label, to file 'roll.dat'
        #{$prg} #{sub} -e 1 -p chi1
             # extract the chi torsion angle of strand I, but exclude
             # those from the two terminal base pairs. For comparison,
             # run also: #{$prg} #{sub} -p chi1
        #{$prg} #{sub} -a
             # extract all parameters, each in a separate file
Options:
#{$sline}
TXT
        end

        def parse_end_bps(end_bps)
            ib, ie = end_bps
            ib = 0 if ib < 0
            ie ||= ib
            ie = 0 if ie < 0
            $stderr.puts "\tend effects set to: -e #{ib} #{ie}" unless [ib, ie] == [0, 0]
            [ib, ie]
        end

        def parse_options
            txt = docinfo
            opts = Trollop::options do
                banner txt
                opt :separator, "Separator for fields [\\t]", :default => "\t"
                opt :par_name, "Name of parameter to extract", :type => :string
                opt :fromfile, "Parameters file", :default => "#{$oname}.out"
                opt :outfile, "File of selected parameter", :default => "stdout"
                opt :end_bps, "Number of end pairs to ignore", :default => [0, 0]

                opt :all, "Extract all parameters into separate files", :default => false
                opt :clean, "Clean up parameter files by the -a option", :default => false
                opt :list, "List all parameters", :default => false
                opt :no_1col, "Delete the first (label) column", :default => false
            end
            opts[:separator] = opts[:separator][0, 1] # keep just the first character
            opts[:end_bps] = parse_end_bps(opts[:end_bps])
            opts[:all] = false if opts[:par_name]

            Trollop::die :fromfile, "must exist" unless File.exist?(opts[:fromfile])

            opts
        end
    end

    class Reorient
        def main
            Ensemble.set_globals
            opts = parse_options

            if opts[:ensemble]
                ensemble = opts[:ensemble]
                models = Utils.extract_model_numbers_from_ensemble_pdb(ensemble)
            elsif opts[:models]
                ensemble, models = Utils.readlist_model_numbers(opts[:models])
            elsif opts[:list]
                pdblist = Utils.collect_pdb_list(opts[:list])
            else    # use a pattern: "pdbdir/*.pdb", should be quoted
                pdblist = Dir[opts[:pattern]]
            end

            if models
                Utils.fatal("no MODEL/ENDMDL pairs in ensemble file") if models.size == 0
                Utils.output_model_list(ensemble, models, $model_list)
                reorient_ensemble_models(ensemble, models, opts)
            else
                Utils.fatal("no matched PDB files") if pdblist.size == 0
                Utils.output_pdb_list(pdblist, $pdb_list)
                reorient_pdb_list(pdblist, opts)
            end
        end

        def reorient_each_model(opts, renew_bpfile=false)
            if opts[:frame]     # for frame_mol
                if opts[:single]
                    Utils.run_cmd("#{$x3dna}find_pair -s #{$temp_model}.pdb #{$temp_model}.inp")
                else
                    Utils.update_bpfile(opts[:bpfile], $temp_model) if renew_bpfile
                end
                Utils.run_cmd("#{$x3dna}analyze #{$temp_model}.inp #{$msg}") # get ref_frames.dat
                opt = "#{opts[:frame]} ref_frames.dat #{$temp_model}.pdb #{$temp_model}_trx.pdb"
                Utils.run_cmd("#{$x3dna}frame_mol #{opt} #{$msg}")
            else                # rotate_mol
                opt = "#{opts[:rotate]} #{$temp_model}.pdb #{$temp_model}_trx.pdb"
                Utils.run_cmd("#{$x3dna}rotate_mol #{opt} #{$msg}")
            end
        end

        def append_transformed_coordinates(aFile, model, pdbfile, hmsg)
            aFile.puts hmsg unless hmsg.empty?
            aFile.puts "%6s    %4d" % ["MODEL ", model]
            aFile.puts File.read(pdbfile).sub(/^END\s*$/, "ENDMDL")
        end

        def reorient_ensemble_models(ensemble, models, opts)
            if opts[:keep]    # get rotmat of the representative model
                Utils.get_best_model(ensemble, models, "best_str.pdb")
                FileUtils.cp("best_str.pdb", "#{$temp_model}.pdb")
                reorient_each_model(opts, true)
                FileUtils.cp("rotmat.dat", "best_rotmat.dat")
            end

            tnum = models.size
            File.open(opts[:outfile], "w") do |aFile|
                models.each_with_index do |model, num|
                    puts "Process model ##{model}\t#{num + 1} / #{tnum}"
                    cmd = "#{$x3dna}ex_str -nmr -#{model} #{ensemble} #{$temp_model}.pdb"
                    Utils.run_cmd("#{cmd} #{$msg}")
                    if opts[:keep]
                        opt = "-t=best_rotmat.dat #{$temp_model}.pdb #{$temp_model}_trx.pdb"
                        Utils.run_cmd("#{$x3dna}rotate_mol #{opt} #{$msg}")
                    else
                        reorient_each_model(opts, num == 0)
                    end
                    append_transformed_coordinates(aFile, model, "#{$temp_model}_trx.pdb", "")
                end
            end
        end

        def reorient_pdb_list(pdblist, opts)
            if opts[:keep]      # get rotmat of the first model
                FileUtils.cp(pdblist.fitst, "#{$temp_model}.pdb")
                reorient_each_model(opts, true)
                FileUtils.cp("rotmat.dat", "best_rotmat.dat")
            end

            tnum = pdblist.size
            File.open(opts[:outfile], "w") do |aFile|
                pdblist.each_with_index do |pdb, num|
                    puts "Process PDB file #{pdb}\t#{num + 1} / #{tnum}"
                    FileUtils.cp(pdb, "#{$temp_model}.pdb")
                    if opts[:keep]
                        opt = "-t=best_rotmat.dat #{$temp_model}.pdb #{$temp_model}_trx.pdb"
                        Utils.run_cmd("#{$x3dna}rotate_mol #{opt} #{$msg}")
                    else
                        reorient_each_model(opts, num == 0)
                    end
                    hmsg = "REMARK ------ based on PDB file '#{pdb}'"
                    append_transformed_coordinates(aFile, num + 1, "#{$temp_model}_trx.pdb", hmsg)
                end
            end
        end

        def docinfo
            sub = 'reorient'
            <<TXT
#{$sline}
Reorient models of a MODEL/ENDMDL delineated ensemble of NMR structures
or MD trajectories based on user specified base-pair reference frame or
rotation matrix. Coordinate transformation is perform by 'frame_mol' or
'rotate_mol'. Useful for 'local' structural alignment. The transformed
ensemble can be converted to an image using '#{$prg} block_image'
and visualized with Jmol or PyMOL.

Usage:
        #{$prg} #{sub} options
Examples:
        #{$prg} #{sub} -f 'm -6' -b bpfile.dat -e 2kei.pdb
             # reorient each model using the reference frame of bp #6
             #   with minor groove facing the viewer (-m)
             # options to -f must be quoted, where dash is optional
             #   generate '#{$model_list}', '#{$oname}_trx.pdb'
        #{$prg} #{sub} -f '8,9' -b bpfile.dat -e 2kei.pdb
             # reorient each model using the middle reference frame
             #   between base-pairs #8 and #9
        #{$prg} #{sub} -v -e 2kei.pdb -o extended_view.pdb
             # transform the ensemble to its most extended view and
             #   keep the original relative orientation.
Options:
#{$sline}
TXT
        end

        def parse_options
            txt = docinfo
            opts = Trollop::options do
                banner txt
                opt :bpfile, "Name of file containing base-pairing info", :type => :string
                opt :outfile, "Output file", :default => "#{$oname}_trx.pdb"

                opt :single, "Single-stranded DNA/RNA", :default => false

                opt :ensemble, "Ensemble delineated with MODEL/ENDMDL pairs", :type => :string
                opt :models, "File containing an explicit list of model numbers", :type => :string
                opt :pattern, "Pattern of model files to process (e.g., *.pdb)", :type => :string
                opt :list, "File containing an explicit list of models", :type => :string

                opt :rotate, "Options to be transfered to rotate_mol (quoted)", :type => :string
                opt :frame, "Options to be transfered to frame_mol (quoted)", :type => :string

                opt :keep, "Keep the original relative orientation", :default => false
                opt :bestview, "Set to the most extended view; keep orientation", :default => false, :short => '-v'
            end
            k = Utils.count_set_options([:ensemble, :models, :pattern, :list], opts)
            Utils.fatal("Specify only one of: ensemble|models|pattern|list") unless k == 1

            Trollop::die :ensemble, "must exist" if opts[:ensemble] && !File.exist?(opts[:ensemble])
            Trollop::die :models, "must exist" if opts[:models] && !File.exist?(opts[:models])
            Trollop::die :list, "must exist" if opts[:list] && !File.exist?(opts[:list])

            k = Utils.count_set_options([:rotate, :frame], opts)
            if k == 0 || opts[:bestview]
                opts[:rotate] = "-cb" # default: set to the best view on DNA/RNA
                opts[:keep] = true    # same as the original Perl script 'nmr_ensemble'
            elsif k == 2              # frame_mol takes precedence
                puts "ignoring -r option '#{opts[:rotate]}', use -f option '#{opts[:frame]}'"
                opts[:rotate] = nil
            end

            if opts[:frame]
                Utils.assert_file_exist(opts[:bpfile], "must specify a bpfile with the -f option") unless opts[:single]
                opts[:frame] = Utils.check_dashes(opts[:frame])
            else                # opts[:rotate]
                opts[:rotate] = Utils.check_dashes(opts[:rotate])
            end

            opts
        end
    end

    class BlockImage
        def main
            bkview = "blocview" # call the Ruby version

            Ensemble.set_globals
            opts = parse_options
            ensemble = ARGV.shift
            Utils.assert_file_exist(ensemble, "must specify a PDB file")

            models = Utils.extract_model_numbers_from_ensemble_pdb(ensemble)
            Utils.fatal("no MODEL/ENDMDL pairs in ensemble file") if models.size == 0

            selfile = "picked_str.pdb"
            bstr_num = Utils.get_best_model(ensemble, models, selfile)
            puts "Process the best representative model <##{bstr_num}>"
            opt = opts[:bkview]
            Utils.run_cmd("#{$x3dna}#{bkview} #{opt} #{selfile} #{$msg}")

            FileUtils.cp("blocview.r3d", opts[:r3dfile])
            if opts[:repr]      # just the representative model
                FileUtils.cp("blocview.png", opts[:imgfile])
            else                # the whole ensemble
                models.each do |m|
                    next if m == bstr_num # already counted
                    Utils.run_cmd("#{$x3dna}ex_str -#{m} -nmr #{ensemble} #{selfile} #{$msg}")
                    puts "Process model <##{m}>"
                    Utils.run_cmd("#{$x3dna}#{bkview} #{opt} #{selfile} #{$msg}")
                    Utils.remove_r3d_header("blocview.r3d", "blocview_noheader.r3d")
                    Utils.append_r3dfile(opts[:r3dfile], "blocview_noheader.r3d")
                end
                Utils.x3dna_r3d2png(opts)
                Utils.output_scale_factor(opts[:r3dfile]) unless opts[:dpi_pymol]
            end
            File.delete(selfile)
        end

        def docinfo
            sub = 'block_image'
            <<TXT
#{$sline}
Ensemble version of 'blocview', but WITHOUT coordinate transformation.
Run '#{$prg} reorient' to set the structures in your preferred
orientation. Raster3D (or PyMOL) and ImageMagick must be installed.

Usage:
        #{$prg} #{sub} options PDBfile
Examples:
        #{$prg} #{sub} 2kei.pdb
             # the whole ensemble, image named 'ensemble.png'
        #{$prg} #{sub} -m -i 2kei-selected.png 2kei.pdb
             # the best representative model
Options:
#{$sline}
TXT
        end

        def parse_options
            txt = docinfo
            opts = Trollop::options do
                banner txt
                opt :imgfile, "name of image file", :default => "ensemble.png"
                opt :r3dfile, "name of .r3d file", :default => "ensemble.r3d"

                opt :dpi_pymol, "create PyMOL ray-traced image at specific DPI", :type => :integer
                opt :scale, "set scale factor (for 'render' of Raster3D)", :type => :float
                opt :repr, "use the representative (or first) model", :default => false, :short => "-m"
                opt :ball_and_stick, "get a ball-and-stick image", :default => false
                opt :p_base_ring, "use only P and base ring atoms", :default => false, :short => "-c"
                opt :no_ds, "do not show double-helix ribbon", :default => false
            end
            opts[:imgfile] += ".png" unless opts[:imgfile] =~ /\.png$/i

            bkview = "--original"
            bkview += " --dpi-pymol #{opts[:dpi_pymol]}" if opts[:dpi_pymol]
            bkview += " --scale #{opts[:scale]}" if opts[:scale]
            bkview += " --ball-and-stick" if opts[:ball_and_stick]
            bkview += " --p-base-ring" if opts[:p_base_ring]
            bkview += " --no-ds" if opts[:no_ds]
            opts[:bkview] = bkview # to be transfered to 'blocview'

            opts
        end
    end
end
