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
#  Module: Parser::SHORE::Variants::select_by_quality.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "\n$0 shore_peak_summary\n\n";
my $file = shift or die $usage;

open FILE, $file or die $usage;

while (<FILE>) {

	chomp;

	if( substr($_, 0, 1) ne "#" ) {

		my @a = split(/\t/, $_);

		my $beg = $a[2] - 1;
		my $end = $beg + $a[3] - 1;

		print $a[0] ."\tchr". $a[1] ."\t$beg\t$end\t+\t$a[4]\t$a[6]\t$a[7]\t$a[8]\t0\t$a[5]\n";
	}
}



