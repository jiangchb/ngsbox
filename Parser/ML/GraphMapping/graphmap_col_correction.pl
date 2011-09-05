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
#  Module: Parser::ML::GraphMapping::graphmap_col_correction.pl
#  Purpose:
#  In:
#  Out:
#



# The transformation of Col-0 alignments was wrong:
# Just copy the Strain Col-0 alignment to the Reference alignment if different

my $usage = "$0 map.list\n";
my $file = shift or die $usage;
open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	my $eco = $a[0];
	if ($eco eq "Col-0") {
		if ($a[4] eq $a[5]) {
			print $line;
		}
		else {
			print $a[0], "\t", $a[1], "\t", $a[2], "\t", $a[3], "\t", $a[5], "\t", $a[5];
			for (my $i = 6; $i < @a+0; $i++) {
				print "\t", $a[$i];
			}
			print "\n";
		}
	}
	else {
		print $line;
	}
}

