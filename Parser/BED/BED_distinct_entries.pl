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
#  Module: Parser::BED::BED_distinct_entries.pl
#  Purpose:
#  In:
#  Out:
#


my $extension = shift;# or die;
my $file = shift or die;
open IN, $file or die "Cannot open input file\n";

my $last_chr = "NA";
my $last_beg = -999999;
my $last_end = -999999;


while( <IN> ) {

	chomp;
	my ($chr, $beg, $end, $name) = split(/\t/, $_);

	if($chr ne $last_chr || $beg != $last_beg || $end != $last_end) { 
		$beg -= $extension;
		$end += $extension;
		if($beg < 1) { $beg = 1; }

		print "$chr\t$beg\t$end\t$name\n";
	}
	
	$last_chr = $chr;
	$last_beg = $beg;
	$last_end = $end;
}

close IN;

exit(0);
