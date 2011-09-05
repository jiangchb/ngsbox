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
#  Module: Parser::FL::prb2sanger.pl
#  Purpose:
#  In:
#  Out:
#


my $file = shift;

my $qual_offset = 0;
my $qual2sanger = 1;

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;
	my ($id, $sequence, $pe, $qual1, $qual2, $qual3) = split(/\t/, $_);

	if ($qual_offset == 1) {
		my $new = "";
		for (my $i = 0; $i<length($qual1); $i++) {
			$new .= chr((ord(substr($qual1, $i, 1))+14));
		}
		$qual1 = $new;
	}

	if ($qual2sanger == 1) {
		$qual2 = "";
		for (my $i = 0; $i<length($qual1); $i++) {
			$qual2 .= chr(int(33 + 10 * log(1 + 10**((ord(substr($qual1, $i, 1)) - 64)/10.0)) / log(10)+0.499 ));
		}
	}

	print "$id\t$sequence\t$pe\t$qual1\t$qual2\t$qual3\n";
}

close IN;

exit(0);
