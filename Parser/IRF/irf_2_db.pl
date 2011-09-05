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
#  Module: Parser::IRF::irf_2_db.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 irfoutput";
my $file = shift or die $usage;

open FILE, $file or die $usage;

my $chr = -1;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (@a > 15) {
		print $chr;
		for (my $i = 0; $i < @a; $i++) {
			print "\t", $a[$i];
		}
		print "\n";
	}
	else {
		if (substr($line, 0, 8) eq "Sequence") {
			$chr = @a[1];
print STDERR "Now parsing chromsome:", $chr, "\n";
		}
	}
}

