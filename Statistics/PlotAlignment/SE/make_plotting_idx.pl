#!/usr/bin/perl
####################################################################################
#Author 	Korbinian Schneeberger 
#Date 		05/15/07
#Version	0.1
#Input		map.list
#Function	Index of a map.list file for fast access
####################################################################################

use strict;
use warnings;
use Getopt::Long;

my $file;
my $scale;
my %CMD;

GetCom();

open FILE, $file or die "Cannot open $file\n";
open OUT, "> $file.idx" or die "Cannot open idx file\n";

my $last_chr = 1;
my $last_tell = 0;
my $last_pos = 0; 
#my $last_print_chr = -1;
#my $last_print_pos = -1;
my $lc = 0;

while(<FILE>) {
	$lc++;
	my @a = split " ";
	if ($lc % $scale == 0) {
		#if ($last_print_chr != $last_chr or $last_print_pos != $last_pos) {
			print OUT $last_chr, "\t", $last_pos, "\t", $last_tell, "\n";
		#}
		#$last_print_chr = $last_chr;
		#$last_print_pos = $last_pos;
	}
	$last_tell = tell(FILE);
	$last_chr = $a[0];
	$last_pos = $a[1];
}


sub GetCom {

  my @usage = ("\nUsage: $0 --file=<file> --scale=<int>

required:
--file\t\tmap.list
--scale\t\tDistance btw marks
\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "file=s", "scale=s");

	die("Please specify an input file\n") unless $CMD{file};
	die("Please specify scale\n") unless $CMD{scale};

	$file = $CMD{file};
	$scale = $CMD{scale};

	return(0);
}


exit(0);
