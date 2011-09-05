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
#  Module: Analysis::Assembly::TargetedAssemblyVelvet::reassembly_WG.pl
#  Purpose:
#  In:
#  Out:
#

###########################################################################
# Author 	Korbinian Schneeberger 
# Date 		09/26/07
# Version	0.2
# Input		map.list, map.list.idx, chr, start, end
# Function	returns all reads from a specified region of the mapping
###########################################################################


use lib "/ebio/abt6/stephan/pgsp/Assembly";
use parse_mapping;

my $file	= shift;
my $chr		= shift;
my $begin	= shift;
my $end		= shift;
my $read_length	= shift;
my $out_dir	= shift;


for(my $i = $begin; $i< $end; $i+=100) {
	my $end = $i + 200;
	my @results = parse_mapping::get($file, $chr, $i - $read_length + 1, $end);

	open OUT, ">$out_dir/$chr-$i-$end.fa" or die;

	foreach (@results) {
		my @entries = split(" ", $_);
		my $seq = $entries[2];
	
		print OUT ">$entries[3] | $entries[0] | $entries[1]\n";

		for(my $i = 0; $i < length($seq); $i++) {
			if(substr($seq, $i, 1) eq "[") {
				$i+=2;
				if(substr($seq, $i, 1) ne "-") {
					print OUT substr($seq, $i, 1);
				}
				$i++;
			}
			else {
				if(substr($seq, $i, 1) ne "-") {
					print OUT substr($seq, $i, 1);
				}
			}
		}
		print OUT "\n";
	}

	close OUT;
}

exit(0);

