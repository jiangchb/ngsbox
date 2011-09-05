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
#  Module: Analysis::Assembly::TargetedAssemblyVelvet::parse_mapping.pm
#  Purpose:
#  In:
#  Out:
#

###############################################################################
# Author 	Korbinian Schneeberger, Stephan Ossowski 
# Date 		09/26/07
# Version	0.3
# Input		map.list, chr, start, end
# Function	returns all reads starting in a specified region of the mapping
###############################################################################


package parse_mapping;
 
sub get
{
	my ($file, $chr, $begin, $end) = @_;
	my @results = ();
	open FILE, $file or die "Cannot open $file\n";
	open IDX, "$file.idx" or die "Cannot open $file.idx\n";

	### Read in index file an jump FILE to right position
	my $jump = 0;
	SETTING: while ( <IDX> ) {
		my @a = split(" ", $_);
		if ($a[0] > $chr or ($a[0] == $chr and $a[1] >= $begin)) {
			last SETTING;
		}
		else { $jump = $a[2]; }
	}
	seek(FILE, $jump, 0);
	
	### Set file to region of interest and read in first fragment
	FIND: while (<FILE>) {
		my ($current_chr, $current_pos) = split (" ", $_);
		if ($current_chr == $chr and $current_pos>= $begin) {
			chomp;
			push @results, $_;
			last FIND;
		}
	}

	### Parse until the end of region of interest
	REG: while (<FILE>) {
		chomp;
		my ($current_chr, $current_pos) = split (" ", $_);
		if ($chr != $current_chr or $current_pos > $end) {
			last REG;
		}
		push @results, $_;
	}

	return(@results);
}
1;
