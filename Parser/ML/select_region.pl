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
#  Module: Parser::ML::select_region.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "\nselect_region.pl map.list chr begin end\n\n";
my $file = shift or die $usage;
my $chr = shift or die $usage;
my $begin = shift or die $usage;
my $end = shift or die $usage;

open FILE, $file or die "Cannot open file\n";
open OUT, ">".$file.".$chr:$begin..$end" or die "Cannot open file\n";

my $written = 0;

while (my $line = <FILE>) {
	my @a = split " ", $line;

	if ($a[1]%200000 == 0) {
		print STDERR $a[0], "\t", $a[1], "\n";
	}

	if ($a[0] == $chr && $a[1] >= $begin && $a[1] <= $end) {
		print OUT $line;
		$written = 1;
	}
	else {
		if ($written == 1) {
			exit(0);
		}
	}
}


