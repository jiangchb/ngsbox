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
#  Module: Parser::FL::cut_length.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 length file\n";
my $length = shift or die $usage;
my $file   = shift or die $usage;

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;
	my @elements = split(/\t/, $_);
	
	if( length($elements[1]) >= $length ) {
		my $sequence = substr($elements[1], 0, $length);
		my $prb      = substr($elements[3], 0, $length);
		my $qCal     = substr($elements[4], 0, $length);
		my $chas     = substr($elements[5], 0, $length);
		
		print "$elements[0]\t$sequence\t$elements[2]\t$prb\t$qCal\t$chas\n";
	}
}

close IN;

exit(0);
