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
#  Module: Support::Misc::get_columns.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 file\n";
my $file  = shift or die $usage;

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;

	my @e = split("\t", $_);

	my $sample = "";
	if(length($e[0]) == 1) {
		$sample = "00" . $e[0] . "TD";
	}
	elsif(length($e[0]) == 2) {
		$sample = "0" . $e[0] . "TD";
	}
	elsif(length($e[0]) == 3) {
		$sample = $e[0] . "TD";
	}

	print "$sample\t$e[1]\t$e[2]\t$e[3]\t$e[4]\n";
}

close IN;

exit(0);
