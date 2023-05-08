#!/usr/bin/env perl
## a simple script for setting up 3DNA environment

use warnings;
use strict;
use FindBin qw($Bin);
use File::Basename;

my $X3DNA_Dir = $Bin;
$X3DNA_Dir =~ s/\/\w+$//;

my $shell    = $ENV{'SHELL'} ? basename( $ENV{'SHELL'} ) : 'bash';

my $num_dash = 78;
my ( $rcfile, $env_mtx, $path_mtxbin );

if ( $shell =~ /bash/ ) {
    $rcfile      = "~/.bashrc";
    $env_mtx     = "export X3DNA='$X3DNA_Dir'";
    $path_mtxbin = "export PATH='$X3DNA_Dir/bin':\$PATH";
} elsif ( $shell =~ /csh/ ) {
    $rcfile      = "~/.tcshrc or ~/.cshrc";
    $env_mtx     = "setenv X3DNA '$X3DNA_Dir'";
    $path_mtxbin = "setenv PATH '$X3DNA_Dir/bin':{\$PATH}";
} else {
    print "Sorry, I am not familiar with your $shell shell!\n";
    $rcfile      = "unknown";
    $env_mtx     = "unknown";
    $path_mtxbin = "unknown";
}

my $setup = <<X3DNA_SETUP;

To install X3DNA, do as follows:
  (0) download 3DNA binary distribution for your system from URL
          http://x3dna.org
  (1) tar pzxvf x3dna_v2.1_Linux.tar.gz
  (2) cd x3dna_v2.1/bin
  (3) ./x3dna_setup
        To run X3DNA, you need to set up the followings:
          o the environment variable X3DNA
          o add \$X3DNA/bin to your command search path
        for your '$shell' shell, please add the following into $rcfile:
          --------------------------------------------------------------
            $env_mtx
            $path_mtxbin
          --------------------------------------------------------------
  (4) type find_pair -h
        for command line help -- this applies to all 3DNA binaries
        Visit 3DNA homepage at URL http://x3dna.org/ for more info.
X3DNA_SETUP

print '+' x $num_dash;
print $setup;
print '+' x $num_dash, "\n";
