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
#  Module: Analysis::CNV::join_regions.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 duplication.txt\n";

open FILE, shift or die $usage;

my $sample = "";
my $chr = -1;
my $start = -1;
my $end = -1;
my $cnv_count = 0;

while (<FILE>) {
	my @a = split " ";
	if ($chr != $a[1] or $a[2] > $end + 5000) {
		if ($chr != -1) {
			print $sample, "\t", $chr, "\t", 	$start, "\t", $end, "\t", $cnv_count, "\n";
		}
		$cnv_count = 0;
		$start = $a[2];
	}

	$sample = $a[0];
	$chr = $a[1];
	$end = $a[3];
	$cnv_count += $a[5];
}

if ($chr != -1) {
	print $sample, "\t", $chr, "\t",  $start, "\t", $end, "\t", $cnv_count, "\n";
}



