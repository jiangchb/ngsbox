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
#  Module: Support::Misc::pos_2_region.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 file(chr pos)\nWill be translated in (chr begin end)\n";
my $file = shift or die $usage;

open FILE, $file or die "Cannot open file\n";

my $bc = -1;
my $bp = -1;
my $lp = -1;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[0] > $bc or $a[1] > $lp+1) {
		if ($bc != -1) {
			print $bc, "\t", $bp, "\t", $lp, "\n";
		}
		$bc = $a[0];
		$bp = $a[1];
	}
	$lp = $a[1];
}
if ($bc != -1) {
	print $bc, "\t", $bp, "\t", $lp, "\n";
}

