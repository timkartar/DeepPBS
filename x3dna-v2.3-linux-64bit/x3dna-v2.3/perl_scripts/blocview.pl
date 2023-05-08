#!/usr/bin/env perl
use warnings;
use strict;
use vars qw($opt_o $opt_d $opt_b $opt_r $opt_p $opt_a $opt_h
    $opt_t $opt_c $opt_C $opt_k $opt_s $opt_i $opt_j $opt_x
    $opt_y $opt_z $X3DNA_BIN $dft_r3dfile $dft_pdbfile);

use Getopt::Std;    # for easy command-line options processing
use File::Copy;     # File copy

use FindBin qw($Bin);
use lib $Bin;
use x3dna_utils;

process_options();

$X3DNA_BIN = x3dna_bin();
system("$X3DNA_BIN/../perl_scripts/dcmnfile -x") unless $opt_C;
die "\n" if ($opt_c);    # just for file clean up

$dft_r3dfile = "t.r3d";
$dft_pdbfile = "t.pdb";

output_usage() if ( @ARGV != 1 or $opt_h );

my $pdbfile = validate_pdbfile( $ARGV[0], $dft_pdbfile );
copy( $pdbfile, $dft_pdbfile );    # make a copy to a default file name
my $imgfile = ($opt_j) ? "t.jpg" : "t.png";

set_view();
get_r3d_header();
get_component_scenes();
combine_scenes();

render_image( $dft_r3dfile, $imgfile, $opt_t );
copy( $imgfile, $opt_i ) if ($opt_i);    # save image to its designated name

output_scale_factor($dft_r3dfile) unless ($opt_t);

# get special PDB file name ($dft_pdbfile) back
copy( $pdbfile, $dft_pdbfile ) if ( $pdbfile eq '__copy__' );

# finally display the image if -d option is given
system( "display " . $imgfile ) if ($opt_d);

### --------------------------------------------------------------------
sub combine_scenes {    # combine header (t0.r3d) and all 5 parts together

    open( FOUT, ">$dft_r3dfile" )
        || die "Cann't open file <$dft_r3dfile> for writing: $!\n";
    for ( my $i = 0 ; $i <= 5 ; $i++ ) {
        next if ( $opt_b && $i >= 1 && $i <= 3 );
        my $nfile = "t$i.r3d";
        if ( -e $nfile ) {
            open( FINP, $nfile ) || die "Cann't open file <$nfile> for reading: $!\n";
            print FOUT while (<FINP>);
            close(FINP);
        }
    }
    close(FOUT);
}

sub get_component_scenes {

    if ($opt_b) {    # ball-and-stick model
        copy( $dft_pdbfile, "temp" );
    } else {
        if ($opt_r) {    # P + base Ring atoms
            system("$X3DNA_BIN/get_part -z $dft_pdbfile tb.pdb");
        } else {   # get the nucleic acid backbone (without O1P/O2P) and base ring atoms
            system("$X3DNA_BIN/get_part -x $dft_pdbfile tb.pdb");
        }

        # nucleic acid base and sugar are colored by residue type
        system("$X3DNA_BIN/r3d_atom -ncz -r=0.16 tb.pdb t1.r3d");

        # block representation of the bases
        system("$X3DNA_BIN/pdb2img -rcnm $dft_pdbfile t2.r3d");

        # use "molauto" to get protein secondary structure & nucleic acid backbone
        system("molauto -nocentre -ss_hb $dft_pdbfile > temp")
            ;      # H-bond based secondary structure
        system("$X3DNA_BIN/get_part -c $dft_pdbfile temp2");

        remove_nt_backbone("temp2") if $opt_k;

        # use "molscript" to get an input to Raster3D, delete its header section
        system("molscript -r < temp2 | $X3DNA_BIN/../perl_scripts/del_ms -n > t3.r3d");

        # get input for ligands
        system("$X3DNA_BIN/get_part -t $dft_pdbfile temp");
    }

    if (has_coordinates("temp")) {
        system("$X3DNA_BIN/r3d_atom -r=0.0 -b=0.20 -o -n temp t4.r3d");    # ball for atoms
        system("$X3DNA_BIN/r3d_atom -r=0.12 -gn temp t5.r3d");             # gray rod
    }
}

sub has_coordinates {
    my $pdbfile     = shift;
    my $atom_hetatm = 0;
    open( FP, $pdbfile ) || die "Cann't open file <$pdbfile> for reading: $!\n";
    while ( my $str = <FP> ) {
        if ( $str =~ /^(ATOM  |HETATM)/ ) {
            $atom_hetatm = 1;
            last;
        }
    }
    close(FP);
    return $atom_hetatm;
}

sub remove_nt_backbone {

    my $msc0 = shift;
    my $tmsc = "temp.msc";

    open( FINP, $msc0 )    || die "Cann't open file <$msc0> for reading: $!\n";
    open( FMSC, ">$tmsc" ) || die "Cann't open file <$tmsc> for writing: $!\n";

    while ( my $str = <FINP> ) {
        if ( $str =~ /^\s*double-helix/ ) {
            print FMSC "! $str";
        } else {
            print FMSC $str;
        }
    }

    close(FINP);
    close(FMSC);

    copy( $tmsc, $msc0 );
}

sub get_r3d_header {    # for the whole structure

    if ($opt_s) {       # user supplied scale
        system("$X3DNA_BIN/r3d_atom -s=$opt_s $dft_pdbfile $dft_r3dfile");
    } else {            # use whole structure for default scale
        system("$X3DNA_BIN/r3d_atom $dft_pdbfile $dft_r3dfile");
    }

    open( FINP, "$dft_r3dfile" )
        || die "Can't open file <$dft_r3dfile> for reading: $!\n";
    open( FOUT, ">t0.r3d" ) || die "Can't open file <t0.r3d> for writing: $!\n";
    while (<FINP>) {
        print FOUT if ( $. <= 20 );
    }
    close(FINP);
    close(FOUT);
}

sub set_view {

    # get only protein, nucleic acid, ligands and delete H-atoms
    system("$X3DNA_BIN/get_part -pnt $dft_pdbfile temp");
    system("$X3DNA_BIN/get_part -d temp $dft_pdbfile");

    return if ($opt_o);    # using original PDB file

    if ($opt_p) {          # protein part
        system("$X3DNA_BIN/rotate_mol -cp $dft_pdbfile temp");
    } elsif ($opt_a) {     # all atoms
        system("$X3DNA_BIN/rotate_mol -ca $dft_pdbfile temp");
    } else {               # nucleic acid part
        system("$X3DNA_BIN/rotate_mol -cb $dft_pdbfile temp");
    }

    # further rotation by x-, y- or z-axis
    my $rotfile = "rotxyz.ang";
    open( FOUT, ">$rotfile" ) || die "Can't open file <$rotfile> for writing: $!\n";
    print FOUT "by rotation x $opt_x\n";
    print FOUT "by rotation y $opt_y\n";
    print FOUT "by rotation z $opt_z\n";
    close(FOUT);

    system("$X3DNA_BIN/rotate_mol -c -r=$rotfile temp $dft_pdbfile");
}

sub output_usage {

    die <<"BLOCVIEW_USAGE";
===========================================================================
SYNOPSIS
    blocview [OPTION]... PDBFILE
DESCRIPTION
    Generates a schematic image which combines base block representation
    with protein ribbon. The image has informative color coding for the
    nucleic acid part and is set in the "best-view" by default. Users need
    to have MolScript, Raster3D and ImageMagick properly installed on their
    system.
        -o   use original coordinates contained in the PDB data file
        -j   output image in JPG format (default to PNG)
        -t[=]RESOLUTION   create PyMOL ray-traced image at RESOLUTION
        -d   display the generated image using "display" of ImageMagick
        -b   ball and stick model with filled base ring
        -c   clean up temporary common files
        -r   only backbone P atoms + base Ring atoms of nucleic acids
        -p   set the best view based on Protein atoms
        -a   set the best view based on All atoms
        -s[=]NUM    set scale factor for the image
        -i[=]IMAGE  set image file name (default to t.png)
        -x|y|z=ANGLE  rotation around x, y, or z-axis by ANGLE degrees
        PDBFILE   a PDB data file name (other than '$dft_pdbfile')
EXAMPLES
    blocview -d -i=bdl084.png bdl084.pdb
AUTHOR
    $author
    $url
===========================================================================
BLOCVIEW_USAGE
}

sub process_options {

    getopts("Cckojbdrpahs:i:x:y:z:t:");
    $opt_h ||= 0;    # help message
    $opt_o ||= 0;    # original coordinates
    $opt_d ||= 0;    # if to display the image
    $opt_j ||= 0;    # if to output image in JPEG format
    $opt_b ||= 0;    # display ball-and-stick model
    $opt_r ||= 0;    # only P atoms + base-Ring atoms
    $opt_p ||= 0;    # Protein atoms
    $opt_a ||= 0;    # All atoms
    $opt_c ||= 0;    # if just to clean up temporary files
    $opt_C ||= 0;    # no clean up with "dcmnfile -x"
    $opt_k ||= 0;    # don't draw DNA/RNA: 'double-helix chain A;'

    $opt_a = 0 if ( $opt_p && $opt_a );    # protein preferred to all atoms

    $opt_s ||= 0;                          # set scale factor: -s=NUM or -s NUM
    remove_equal_sign($opt_s);
    $opt_s = 0 if ( $opt_s < 0 );

    $opt_i ||= 0;                          # output image name
    remove_equal_sign($opt_i);

    $opt_x ||= 0;                          # Rotation axis and angle
    remove_equal_sign($opt_x);

    $opt_y ||= 0;
    remove_equal_sign($opt_y);

    $opt_z ||= 0;
    remove_equal_sign($opt_z);

    $opt_t ||= 0;                          # create PyMOL ray_traced image
    remove_equal_sign($opt_t);

    check_image_options( $opt_i, $opt_j );
}
