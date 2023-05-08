module FileList
    CMN_CFILES = %w{
            nrutil
            reb_fncs
            ana_fncs
            app_fncs
            cmn_fncs
            fncs_slre
        }.sort
    
    PRG_CFILES = %w{
            alc2img
            cehs
            fiber
            stack2img
            find_pair
            o1p_o2p
            rebuild
            std_base
            analyze
            comb_str
            frame_mol
            pdb2img
            regular_dna
            step_hel
            anyhelix
            ex_str
            get_part
            r3d_atom
            rotate_mol
            mutate_bases
            find_platform
            lsfit_2strs
            snap
            get_radius
            find_poco
            pdbml2pdb
            backbone
        }.sort

    LUX_PRIVATE = %w{
            lsfit_2strs
            snap
            get_radius
            find_poco
            pdbml2pdb
            backbone
            find_platform
        }.sort

    DIST_PRGS = PRG_CFILES - LUX_PRIVATE

    ALL_CFILES = (CMN_CFILES + PRG_CFILES).sort
    
    HEADER_FILES = %w{
            x3dna
            x3dna_fncs
        }.sort
end
