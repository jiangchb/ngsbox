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
#  Module: Analysis::Eighties::CNV::parse_4_nbslrr_plot.pl
#  Purpose:
#  In:
#  Out:
#

my $file = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/Col-0-Geneva/AlignmentFolder/CNV_per_segment/shore_count.ALL_80.Bak-7.notcomplete.txt";
open FILE, $file or die "Cannot open file:".$file."\n";
open BACK, ">/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/Col-0-Geneva/AlignmentFolder/CNV_per_segment/shore_count.ALL_80.Bak-7.notcomplete.frac_back.txt";
open NBS, ">/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/Col-0-Geneva/AlignmentFolder/CNV_per_segment/shore_count.ALL_80.Bak-7.notcomplete.frac_nbs.txt";
while (my $line = <FILE>) {
	my @a = split " ", $line;
	for (my $i = 7; $i <= 246; $i+=3) {
		if ($i != 13) { # excl bak-7
		my $count = $a[$i];
		my $rpm = $a[$i+1];
		my $frac = $a[$i+2];
		if ($a[1] eq "NBS_LRR_active") {
			print NBS $frac, "\n"; 
		}
		print BACK $frac, "\n"; 
		}
	}	
}


