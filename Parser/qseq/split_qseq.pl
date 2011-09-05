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
#  Module: Parser::qseq::split_qseq.pl
#  Purpose:
#  In:
#  Out:
#



my $file   = shift or die "$0 file";
my $current_lane = -1;

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;
	my @e = split(/\t/, $_);

	if($current_lane != $e[2]) {
		if($current_lane != -1) { close OUT; }
		$current_lane = $e[2];
		open OUT, ">$file.$current_lane" or die "Cannot open output file";
	}
	print OUT "$_\n";
}

close IN;
close OUT;

exit(0);
