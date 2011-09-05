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
#  Module: Support::Misc::concat_missing_pair.pl
#  Purpose:
#  In:
#  Out:
#


my $win_size = shift;
my $conc_size = shift;
my $file = shift;
open FILE, $file or die "Cannot open infile\n";

my $current_chr = -1;
my $start       = -1;
my $end         = -1;

while(<FILE>) {
	my ($chr, $pos, $count, $freq) = split " ";

	if( ($current_chr != $chr) || ($pos > $start + $conc_size) ) {
		my $length = $end - $start + 1;
		print "$current_chr\t$start\t$end\t$length\n";
		$current_chr = $chr;
		$start = $pos;
	}
	$end = $pos + $win_size;
}

exit(0);
