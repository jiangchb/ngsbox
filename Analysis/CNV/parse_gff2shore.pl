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
#  Module: Analysis::CNV::parse_gff2shore.pl
#  Purpose:
#  In:
#  Out:
#

my $usage = "$0 gff [gene|transposable_element_gene]\n";

my $file = shift or die $usage;
my $feature = shift or die $usage;
open FILE, $file or die "Cannot open file\n";

print "#?chr\tpos\tsize\tstrand\tid\n";

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[2] eq $feature) {
		my @b = split "=", $a[8];
		my $id = substr($b[1], 0, 9);
		my $chr = substr($a[0], 3, 1);
		if ($chr eq "1" or $chr eq "2" or $chr eq "3" or $chr eq "4" or $chr eq "5") {
			print $chr, "\t", $a[3], "\t", $a[4]-$a[3]+1, "\t.\t", $id, "\n";
		}
	}
}


