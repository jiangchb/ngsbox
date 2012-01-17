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
#  Module: Parser::SHORE::Variants::shorePeak2BED.pl
#  Purpose: Convert 'shore peak' output to BED format
#  In: SUMMARY.txt (from shore peak)
#  Out: std::out in BED format
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

		# chr, start, end, name, rank, fake-strand
		# print "chr". $a[1] ."\t$beg\t$end\tchr-$a[1]-$beg\t$a[4]\t+\n";
		
		# chr, start, end, name
		print "chr". $a[1] ."\t$beg\t$end\tpeak_" . $a[4] . "\n";
	}
}



