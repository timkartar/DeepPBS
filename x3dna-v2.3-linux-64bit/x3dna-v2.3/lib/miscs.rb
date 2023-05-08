### Miscellaneous utility functions to replace their Perl counterparts
###     3DNA v2.3-2017feb08, created and maintained by Xiang-Jun Lu (PhD)

$VERBOSE = true
require 'fileutils'

abort "Please set the 3DNA environment variable" unless ENV['X3DNA']

module Utils
    def self.fatal(msg)
        puts msg unless msg.empty?
        exit 1
    end

    def self.x3dna_homedir
        homedir = ENV['X3DNA']
        unless homedir                  # X3DNA environment variable not set
            home = ENV['HOME'] || ENV['HOMEDRIVE']
            if home && Dir.exists?("#{home}/x3dna-v2.3")
                homedir = "#{home}/x3dna-v2.3"
            else
                abort "Please set the X3DNA environment variable."
            end
        end
        (homedir =~ /\/$/) ? homedir : "#{homedir}/" # with ending slash (/)
    end

    def self.x3dna_bindir
        bdir = x3dna_homedir + 'bin/'  # with ending slash
        bdir.include?(' ') ? "" : bdir
    end

    def self.version
        File.read("#{x3dna_homedir}config/version")
    end

    def self.check_duplicates(items)
        raise "not an array" unless items.kind_of?(Array)
        d = items.size - items.uniq.size
        fatal "With #{d} repeated entries" unless d == 0
    end

    def self.run_cmd(cmd, echo=false)
        puts cmd if echo
        $errcmd = false

        unless system(cmd)      # unless %x[#{cmd}].empty?
            puts "Error in running command: '#{cmd}'"
            $errcmd = true
            return if $errlog

            if cmd.include?('>')
                fatal("Rerun the command w/o redirection for details")
            else
                fatal("")
            end
        end
    end

    def self.msgfile
        "2> msgfile"
    end

    # cross-platform way of finding an executable in the $PATH.
    #     based on http://stackoverflow.com/questions/2108727
    def self.has_executable?(cmd)
        exts = ENV['PATHEXT'] ? ENV['PATHEXT'].split(';') : ['']
        ENV['PATH'].split(File::PATH_SEPARATOR).each do |path|
            exts.each do |ext|
                exe = "#{path}/#{cmd}#{ext}"
                return exe if File.executable? exe
            end
        end
        return nil
    end

    def self.skip_lines(aFile, num=1)
        num.times { aFile.gets }
    end

    def self.comment_empty_line?(line)
        line =~ /^\s*(#|$)/
    end

    def self.separation_line(c='-', num=72)
        c * num
    end

    def self.check_dashes(opt)
        # make dash (-) within in 'opt' optional. To transfer options
        # from 'x3dna_ensemble reorient' to either 'frame_mol' or
        # 'rotate_mol', the first dash needs to be escaped: e.g.,
        # -f '\-m -6' , which is easy to get wrong. With the function,
        # it can now be simplified: -f 'm 6' or -f 'm -6'
        # note, with double dash (--), the followings are also fine:
        #     --frame '-m -6'     OR     --frame='-m -6'
        opt.split.map { |o| (o =~ /^-/) ? o : "-#{o}" }.join(" ")
    end

    def self.count_set_options(syms, opts)
        syms.inject(0) { |sum, ele| sum += opts[ele] ? 1 : 0 }
    end

    def self.assert_file_exist(filename, msg="")
        fatal("#{msg}: filename nil") unless filename
        fatal("File '#{filename}' does not exist") unless File.exist?(filename)
    end

    def self.extract_model_numbers_from_ensemble_pdb(ensemble)
        models = []
        n1 = n2 = 0
        File.open(ensemble).each_line do |line|
            if line =~ /^MODEL\s+(\d+)/
                models << $1.to_i
                n1 += 1
            elsif line =~ /^ENDMDL/
                n2 += 1
            end
            fatal("File #{ensemble}: unbalanced MODEL/ENDMDL pairs") unless n1 == n2 || n1 == n2 + 1
        end
        fatal("File #{ensemble}: unmatched MODEL/ENDMDL pairs") unless n1 == n2

        $stderr.puts "\t[i] #{ensemble}: with model numbers <= 0" if models.any? { |x| x <= 0 }
        check_duplicates(models)

        models
    end

    def self.get_best_model(ensemble, models, bestfile)
        bstr_num = models.first # default to the first model
        run_cmd("#{x3dna_bindir}ex_str -nmrb #{ensemble} #{bestfile} #{msgfile}")
        File.open("msgfile").each_line do |line|
            if line =~ /Model selected \#(\d+)/
                bstr_num = $1.to_i
            end
        end
        bstr_num
    end

    def self.readlist_model_numbers(listfile)
        models = []
        ensemble = ""

        File.open(listfile) do |aFile|
            ensemble = aFile.gets.split.first # must be in the first line
            fatal("'#{listfile}': non-existent <#{ensemble}>") unless File.exist?(ensemble)
            while line = aFile.gets do
                next if comment_empty_line?(line)
                models << line.to_i
            end
        end
        check_duplicates(models)

        models_full = extract_model_numbers_from_ensemble_pdb(ensemble)
        dd = models - models_full
        fatal("'#{listfile}' contains non-existent models: <#{dd.join(', ')}>") unless dd.empty?

        [ensemble, models]
    end

    def self.output_model_list(ensemble, models, listfile)
        File.open(listfile, "w") do |aFile|
            aFile.puts ensemble
            aFile.puts <<TXT

# The ensemble PDB file MUST be on the FIRST line. Each model number
# is on a separate line. Empty lines or lines starting with '#' are
# ignored, allowing for easy skip of a model by putting '#' in front
# of it: e.g., #6 to ignore model #6

TXT
            models.each do |model|
                aFile.puts model
            end
        end
    end

    def self.read_pair_info(bpfile)
        num_bp = 0
        bp_info = []
        File.open(bpfile) do |aFile|
            skip_lines(aFile, 2) # excluding the first two lines (.pdb/.out)

            line = aFile.gets
            fatal("Not a duplex: <#{$1}>") unless line =~ /^\s*(\d+)/ && $1.to_i == 2
            bp_info << line

            line = aFile.gets
            if line =~ /^\s*(\d+)/
                num_bp = $1.to_i
                bp_info << line
            else
                fatal("Invalid base-pair file <#{bpfile}>: #{$1} pairs?")
            end

            (1 + num_bp).times { bp_info << aFile.gets }
        end

        [num_bp, bp_info]
    end

    def self.check_bpinfo(bp_info, pdbfile)
        ntsfile = "nts_list.txt"
        run_cmd("#{x3dna_bindir}find_pair -s #{pdbfile} #{ntsfile} #{msgfile}")

        full_nts = []
        File.open(ntsfile) do |aFile|
            skip_lines(aFile, 3)
            num_nts = aFile.gets.to_i
            skip_lines(aFile)   # skip the fifth line
            num_nts.times { full_nts << aFile.gets.to_i }
        end

        bp_nts = []
        bp_info[3..-1].each do |line| # only residue numbers in base-pairs
            bp_nts << line.split[0, 2] # the first two items
        end
        bp_nts = bp_nts.flatten.map { |x| x.to_i }

        dd = bp_nts - full_nts
        fatal("mismatch between base-pair info and model file") unless dd.empty?
    end

    def self.write_pair_info(bp_info, bname)
        File.open("#{bname}.inp", "w") do |aFile|
            aFile.puts "#{bname}.pdb"
            aFile.puts "#{bname}.out"
            bp_info.each { |line| aFile.puts line }
        end
    end

    def self.update_bpfile(bpfile, bname)
        num_bp, bp_info = read_pair_info(bpfile)
        check_bpinfo(bp_info, "#{bname}.pdb")
        write_pair_info(bp_info, bname)

        num_bp
    end

    def self.update_ntlist(bname)
        run_cmd("find_pair -s #{bname}.pdb #{bname}.nts #{msgfile}")
        num_nt = 0
        File.open("#{bname}.nts") do |aFile|
            skip_lines(aFile, 3)
            num_nt = aFile.gets.to_i
        end
        num_nt
    end

    def self.create_pml(pmlfile, r3d_file, pngfile, dpi_pymol)
        File.open(pmlfile, "w") do |aFile|
            aFile.puts "delete all"
            aFile.puts "load #{r3d_file}"

            pymol_file = "pymol_ray.par" # current working directory
            pymol_file = "#{x3dna_homedir}config/#{pymol_file}" unless File.exist?(pymol_file)
            File.open(pymol_file).each_line do |line|
                next if comment_empty_line?(line)
                aFile.puts line
            end

            xpix = dpi_pymol * 8
            aFile.puts "ray #{xpix}"
            aFile.puts "png #{pngfile}"
            aFile.puts "quit"
        end
    end

    def self.collect_pdb_list(listfile)
        pdblist = []
        File.open(listfile).each do |pdb|
            next if comment_empty_line?(pdb)
            pdb.strip!
            if File.exist?(pdb)
                pdblist << pdb
            else
                $stderr.puts "\t[info] ignore -- PDB file '#{pdb}' non-existent"
            end
        end
        check_duplicates(pdblist)

        pdblist                 # keep original order
    end

    def self.output_pdb_list(pdblist, listfile)
        File.open(listfile, "w") do |aFile|
            pdblist.each do |pdb|
                aFile.puts pdb
            end
        end
    end

    def self.remove_r3d_header(inp_r3d, out_r3d)
        File.open(out_r3d, "w") do |oFile|
            File.open(inp_r3d) do |iFile|
                while line = iFile.gets
                    oFile.puts line if $. > 20
                end
            end
        end
    end

    def self.append_r3dfile(dst_r3dfile, src_r3dfile)
        File.open(dst_r3dfile, "a") do |aFile|
            File.open(src_r3dfile).each_line do |line|
                aFile.puts line
            end
        end
    end

    def self.output_scale_factor(r3dfile)
        File.open(r3dfile) do |aFile|
            skip_lines(aFile, 15)
            scale = aFile.gets.split.last.to_f # on the sixth line
            $stderr.puts "***** scale factor used: ===> %.2f <===\n" % scale
        end
    end

    # ****************************** Utilities ******************************
    def self.block_atom_banner
        sub = 'block_atom'
        <<TXT
#{$sline}
Create an ALCHEMY file which combines both atomic and base block
representations, starting from a PDB input file. The output file
has the same basename, but with the ".alc" extension.

Usage:
        #{$prg} #{sub} [-p] pdbfile
Examples:
        #{$prg} #{sub} 355d.pdb
        jmol 355d.alc
        rasmol -alchemy -noconnect 355d.alc
Options:
#{$sline}
TXT
    end

    def self.block_atom(pdbfile)
        fatal("File #{pdbfile} non-existent") unless File.exist?(pdbfile)
        bname = File.basename(pdbfile, '.*')
        bdir = x3dna_bindir

        run_cmd <<CMD
#{bdir}r3d_atom #{pdbfile} temp #{msgfile}
#{bdir}pdb2img #{pdbfile} temp2 #{msgfile}
#{bdir}comb_str atom_lkg.alc bblk_lkg.alc #{bname}.alc #{msgfile}
CMD
    end

    def self.cp_std_banner
        sub = 'cp_std'
        <<TXT
#{$sline}
Select the standard data files to be used with "analyze" and "rebuild".
Available sets include BDNA, ADNA, NDB96 and RNA, which have exactly
the same base geometry and orientation (in the standard base reference
frame) but different backbone conformations. Users can also set their
own standard file with the "std_base" utility program.

Usage:
        #{$prg} #{sub} [options]
Examples:
        #{$prg} #{sub}
             # default to B-DNA backbone conformation
        #{$prg} #{sub} -d ADNA
             # with A-DNA backbone conformation
Options:
#{$sline}
TXT
    end

    def self.cp_std(dataset, block)
        prefix = "#{x3dna_homedir}config"
        if block
            %w(BP R Y).each do |blk| # block with thickness of 1 Angstrom
                FileUtils.cp("#{prefix}/block/Block_#{blk}1.alc", "Block_#{blk}.alc")
            end
        else
            selset = case dataset
                     when /^B/i; "BDNA"
                     when /^A/i; "ADNA"
                     when /^N/i; "NDB96"
                     when /^R/i; "RNA"
                     else
                         puts "Unrecognized dataset: '#{dataset}' -- set to BDNA"
                         "BDNA" # as default"
                     end

            Dir["Atomic[._]?.pdb"].each { |file| File.delete(file) }
            std_bases = Dir["#{prefix}/atomic/#{selset}_?.pdb"] # only standard form
            std_bases.each do |file|
                base = file[-5, 1]
                FileUtils.cp(file, "Atomic_#{base}.pdb")
                lcname = "Atomic.#{base.downcase}.pdb" # lower-case for modified form
                FileUtils.cp(file, lcname)
            end
            FileUtils.cp("#{prefix}/Atomic.p.pdb", "Atomic.p.pdb")
            FileUtils.cp("#{prefix}/Atomic_P.pdb", "Atomic_P.pdb")

            File.open("Atomic_I.pdb", "w") do |aFile|
                num = 0
                File.open("Atomic_G.pdb").each_line do |line|
                    if line =~ /^(ATOM  |HETATM)/
                        next if line =~ / N2 / # skip N2 atom
                        num += 1
                        line.sub!(/  G /, "  I ")
                        line[6, 5] = "%5d" % num
                    end
                    aFile.puts line
                end
            end
            FileUtils.cp("Atomic_I.pdb", "Atomic.i.pdb")
        end
    end

    def self.dcmnfile
        cmnfile = "common_files.dat" # current working directory
        cmnfile = "#{x3dna_homedir}config/#{cmnfile}" unless File.exist?(cmnfile)
        File.open(cmnfile).each_line do |str|
            next if comment_empty_line?(str)
            str.strip!
            File.delete(str) if File.exist?(str)
        end
    end

    def self.x3dna_r3d2png_banner
        sub = 'x3dna_r3d2png'
        <<TXT
#{$sline}
Convert various 3DNA-generated .r3d files to png images using either
Raster3D 'render' or PyMOL ray-tracer. Raster3D/PyMOL and ImageMagick
must be installed.

Usage:
        #{$prg} #{sub} options
Examples:
        find_pair 355d.pdb stdout | analyze
        ex_str -1 bestpairs.pdb bp1.pdb
        r3d_atom -do -r=0.06 -b=0.15 bp1.pdb bp1.r3d
        #{$prg} #{sub} -r bp1.r3d -i bp1_render.png
        #{$prg} #{sub} -d 96 bp1.r3d bp1_pyray.png

        ex_str -5 stacking.pdb bs5.pdb
        stack2img -rtdoc bs5.pdb bs5.r3d
        #{$prg} #{sub} -r bs5.r3d -i bs5_render.png
        #{$prg} #{sub} -d 300 -r bs5.r3d -i bs5_pyray.png
Options:
#{$sline}
TXT
    end

    def self.x3dna_r3d2png(opts)
        convert = "convert"
        opt = "-trim +repage -border 5 -bordercolor white"

        if opts[:dpi_pymol]
            pymol = "pymol"
            if !has_executable?(pymol)
                puts "cannot find '#{pymol}'"
                return false
            end
            cmn = "x3dna_r3d_pymol"
            pmlfile = "#{cmn}.pml"
            tmpimg = "#{cmn}.png"
            create_pml(pmlfile, opts[:r3dfile], tmpimg, opts[:dpi_pymol])
            run_cmd("#{pymol} -qc #{pmlfile} > #{cmn}.msg #{msgfile}")
            if has_executable?(convert)
                run_cmd("#{convert} #{opt} #{tmpimg} #{opts[:imgfile]} #{msgfile}")
            else
                puts "ImageMagick 'convert' cannot be found"
                FileUtils.cp(tmpimg, opts[:imgfile])
            end

        else                    # Raster3D
            render = "render"
            if !has_executable?(render)
                puts "cannot find '#{render}'"
                return false
            end
            tmpimg = "x3dna_r3d.avs"
            run_cmd("#{render} < #{opts[:r3dfile]} > #{tmpimg} #{msgfile}")
            if has_executable?(convert)
                run_cmd("#{convert} #{opt} #{tmpimg} #{opts[:imgfile]} #{msgfile}")
            else
                avsfile = File.basename(opts[:imgfile], ".png")
                puts "ImageMagick 'convert' cannot be found"
                puts "Output image file '#{opts[:imgfile]}' reset to '#{avsfile}'"
                FileUtils.cp(tmpimg, avsfile)
            end
        end

        true                  # with image generated from render/pymol
    end
end
