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
#  Module: Parser::ML::add_pe_column.pl
#  Purpose:
#  In:
#  Out:
#


my $file    = shift;
open IN, $file or die "Cannot open input file\n";
while( <IN> ) {
	chomp;
	my @elements = split(/\t/, $_);
	if(@elements == 12) {
		print "$elements[0]\t$elements[1]\t$elements[2]\t$elements[3]\t$elements[4]\t$elements[5]\t$elements[6]\t$elements[7]\t$elements[8]\t0\t$elements[9]\t$elements[10]\t$elements[11]\n";
	}
	else {
		print "$_\n";
	}
}

close IN;

exit(0);
