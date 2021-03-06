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
#  Module: Analysis::Eighties::SV::Plot::parse_inv_4_plotting.pl
#  Purpose:
#  In:
#  Out:
#

my $usage = "$0 inverionsprediction\n";
my $file = shift or die $usage;
open FILE, $file or die $usage;
my $c = 0;
my $e = "";
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[0] ne $e) {
		$c++;
		$e = $a[0];
	}
	print $c, "\t", $e, "\t", $a[2], "\t", $a[4], "\t", $a[5], "\t", $a[5]-$a[4]+1, "\n"; 
}
