package x3dna_utils;

our @EXPORT = qw(is_integer trim ltrim rtrim x3dna_dir x3dna_bin check_image_options
    remove_equal_sign validate_pdbfile output_scale_factor render_image create_empty_file
    $author $url);
use Exporter;
our @ISA = qw(Exporter);

$author = "3DNA v2.0 [June 8, 2008] (by Dr. Xiang-Jun Lu; 3dna.lu\@gmail.com)";
$url    = "Check http://x3dna.org/ for the latest; Consider upgrade to v2.1";

use strict;
use warnings;

sub is_integer {
    $_[0] =~ /^[+-]?\d+$/;
}

sub trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub ltrim {
    my $string = shift;
    $string =~ s/^\s+//;
    return $string;
}

sub rtrim {
    my $string = shift;
    $string =~ s/\s+$//;
    return $string;
}

sub x3dna_dir {
    my $bdir = $ENV{X3DNA};
    return ( defined $bdir ) ? $bdir : "$ENV{HOME}/X3DNA";
}

sub x3dna_bin {
    return x3dna_dir() . "/bin";
}

sub check_image_options {
    my ( $opt_i, $opt_j ) = @_;

    return unless ( $opt_i and $opt_j );

    die "wrong image type <$opt_i>: must end with '.png' or '.jpg'\n"
        unless ( $opt_i =~ /\.png$/ or $opt_i =~ /\.jpg$/ );

    die "option -j (for jpg image) does not match image name <$opt_i>\n"
        unless ( $opt_i =~ /\.jpg$/ );
}

sub remove_equal_sign {
    $_[0] =~ s/^=//;
}

sub validate_pdbfile {

    my ( $pdbfile, $tfile ) = @_;

    -e $pdbfile or die "PDB file <$pdbfile> does NOT exist!\n";
    if ( $pdbfile eq $tfile ) {    # special case: make a copy
        copy( $pdbfile, "__copy__" );
        $pdbfile = "__copy__";
    }

    return $pdbfile;
}

sub output_scale_factor {
    my $r3dfile = shift;

    open( FINP, $r3dfile ) || die "Can't open file <$r3dfile> for reading: $!\n";
    my $str;
    for ( my $i = 1 ; $i <= 16 ; $i++ ) {
        $str = <FINP>;
    }
    close(FINP);

    my @four = split( /\s+/, trim($str) );
    printf STDERR "***** Scale factor used: ===> %.2f <===\n", $four[3];
}

sub render_image {
    my ( $r3dfile, $imgfile, $pymol_res ) = @_;
    my $pl_dir = x3dna_dir() . "/perl_scripts";

    if ($pymol_res) {    # ray-traced with PyMOL
        qx($pl_dir/x3dna_r3d2png -p -r=$pymol_res $r3dfile $imgfile);
    } else {                            # rendered with Raster3D
        qx($pl_dir/x3dna_r3d2png $r3dfile $imgfile);
    }
}

sub create_empty_file {
    my $filename = shift;

    open( FH, ">$filename" ) or die "Can't open $filename for writing: $!\n";
    close(FH);
}

1;
