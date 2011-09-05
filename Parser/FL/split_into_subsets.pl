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
#  Module: Parser::FL::split_into_subsets.pl
#  Purpose:
#  In:
#  Out:
#



my $split_num = shift;
my $file = shift;

my $count = 0;
my $file_count = 1;

open IN, $file or die "Cannot open input file\n";
open OUT, ">$file_count.txt" or die "Cannot open output file\n";

while( <IN> ) {
	print OUT "$_";
	$count++;

	if($count > $split_num) {
		$count = 0;
		$file_count++;
		close OUT;
		open OUT, ">$file_count.txt" or die "Cannot open output file\n";
	}
}
close IN; close OUT;

exit(0);
