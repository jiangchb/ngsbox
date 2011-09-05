#! /usr/bin/perl
use strict;
use warnings;

###### 
# NGSbox - bioinformatics analysis tools for next generation sequencing data
#
# Copyright 2007-2011 Stephan Ossowski, Korbinian Schneeberger
# 
# NGSbox is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or any later version.
#
# NGSbox is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# Please find the GNU General Public License at <http://www.gnu.org/licenses/>.
#
#  -------------------------------------------------------------------------
#
#  Module: Statistics::PlotAlignment::SE::make_plotting_idx.pl
#  Purpose:
#  In:
#  Out:
#

####################################################################################
#Author 	Korbinian Schneeberger 
#Date 		05/15/07
#Version	0.1
#Input		map.list
#Function	Index of a map.list file for fast access
####################################################################################

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
