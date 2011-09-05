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
#  Module: Parser::ML::prb2sanger.pl
#  Purpose:
#  In:
#  Out:
#


my $file = shift;
open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;
	my @elem = split(/\t/, $_);
	$elem[11] = "";

	for (my $i = 0; $i < length($elem[10]); $i++) {
		$elem[11] .= chr(int(33 + 10 * log(1 + 10**((ord(substr($elem[10], $i, 1)) - 64)/10.0)) / log(10)+0.499 ));
	}

	print "$elem[0]\t$elem[1]\t$elem[2]\t$elem[3]\t$elem[4]\t$elem[5]\t$elem[6]\t$elem[7]\t$elem[8]\t$elem[9]\t$elem[10]\t$elem[11]\t$elem[12]\n";
}

close IN;

exit(0);
