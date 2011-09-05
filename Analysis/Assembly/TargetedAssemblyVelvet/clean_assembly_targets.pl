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
#  Module: Analysis::Assembly::TargetedAssemblyVelvet::clean_assembly_targets.pl
#  Purpose:
#  In:
#  Out:
#


my $flanking_region_length = shift;
my $region_file  = shift;

open IN, "$region_file" or die "Cannot open pc region file\n";

my $last_chr = -1;
my $last_beg = -10000;
my $last_end = -10000;

while(<IN>) {
	chomp;
	my ($chr, $beg, $end) = split (/\t/, $_);

	if( ( $beg > ($last_end + $flanking_region_length) ) || ( $chr != $last_chr) ) {

		my $region_beg = $last_beg - $flanking_region_length;
		my $region_end = $last_end + $flanking_region_length;
		if( $region_beg < 0 ) { $region_beg = 0; }

		if( ( ($region_end - $region_beg < 50000) ) && ($last_chr != -1) ) {
			print "$last_chr\t$region_beg\t$region_end\n";
		}
		$last_beg = $beg;
	}

	if( ($end >= $last_end) || ($chr != $last_chr) ) {
		$last_end = $end;
	}

	$last_chr = $chr;
}

my $region_beg = $last_beg - $flanking_region_length;
my $region_end = $last_end + $flanking_region_length;
if( $region_beg < 0 ) { $region_beg = 0; }

if($region_end - $region_beg < 50000) {
	print "$last_chr\t$region_beg\t$region_end\n";
}

exit(0);
