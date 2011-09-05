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
#  Module: Parser::ML::old_2_new.pl
#  Purpose:
#  In:
#  Out:
#


my $file = shift;

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;
	my @a = split(/\t/, $_);

	my $len = length($a[9]);
	my $qual = "";
	
	for (my $i = 0; $i < $len; $i++) {
		$qual .= "I";
	}

	print $a[0], "\t", $a[1], "\t",$a[2], "\t",$a[3], "\t",$a[4], "\t",$a[5], "\t",$a[6], "\t",$a[7], "\t",$a[8], "\t0\t",$a[9], "\t", $qual, "\t", $a[11], "\n";


}

close IN;

exit(0);
