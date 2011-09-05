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
#  Module: Parser::FL::separate_read_1_and_2.pl
#  Purpose:
#  In:
#  Out:
#



my $file   = shift;

open IN, $file or die "Cannot open input file\n";
open OUT1, ">$file.r1" or die "Cannot open output file\n";
open OUT2, ">$file.r2" or die "Cannot open output file\n";

while( <IN> ) {
	chomp;
	my @e = split(/\t/, $_);
	
	if($e[2] == 1) {
		print OUT1 "$_\n";
	}
	elsif($e[2] == 2) {
		print OUT2 "$_\n";
	}
}

close IN;

exit(0);
