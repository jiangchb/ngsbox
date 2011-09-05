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
#  Module: Parser::ML::gm_score2mm.pl
#  Purpose:
#  In:
#  Out:
#


my $file    = shift;

open IN, $file or die "Cannot open input file\n";

while ( <IN> ) {
	chomp;
	my @e = split(/\t/, $_);

	if(@e != 16) {
		print STDERR "Not enough columns: $e[0]\t$e[1]\n";
	}

	my $brackets = 0;
	for(my $i = 0; $i < length($e[4]); $i++) {

		my $base = substr($e[4], $i, 1);
		if( $base eq "(" || $base eq "[" ) {
			$brackets++;
		}
	}

	print "$e[0]\t$e[1]\t$e[2]\t$e[3]\t$e[4]\t$e[5]\t$e[6]\t$e[7]\t$brackets\t$e[9]\t$e[10]\t$e[11]\t$e[12]\t$e[13]\t$e[14]\t$e[15]\n";
}

close IN;

exit(0);
