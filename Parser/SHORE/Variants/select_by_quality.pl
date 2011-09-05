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


my $usage = "\n$0 quality_variants.txt qthreshold repetitive_threshold support frequency_lb frequency_ub snpsonly\n\n";
my $file = shift or die $usage;
my $qthres = shift or die $usage;
my $rthres = shift or die $usage;
my $support = shift;
my $frequency_lb = shift;
my $frequency_ub = shift;
my $snp = shift;

open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (	($a[8] <= $rthres) and 
		($a[5] >= $qthres) and 
		($snp != 1 || $a[4] ne "-") and 
		($a[7] >= $frequency_lb) and
		($a[7] <= $frequency_ub) and
		($a[6] >= $support)
	) {
		print $line;
	}
}



