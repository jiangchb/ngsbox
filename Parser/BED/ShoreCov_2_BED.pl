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
#  Module: Parser::BED::ShoreCov_2_BED.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "\n$0 min_size segment_file\n\n";
my $min_size = shift or die $usage;
my $segment = shift or die $usage;

open SEG, $segment or die "Cannot open input file\n";

while( my $line = <SEG> ) {
	chomp($line);

	if($line =~ /^#/) {
		#print "$line\n";
	}
	else {
		my @a = split("\t", $line);
		my $beg = $a[1] - 1;
		my $end = $beg + $a[2];
		if($a[2] >= $min_size) {
			print "$a[0]\t$beg\t$end\n";
		}
	}
}

close SEG;

exit(0);
