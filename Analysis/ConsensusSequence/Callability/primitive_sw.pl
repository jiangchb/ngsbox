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
#  Module: Analysis::ConsensusSequence::Callability::primitive_sw.pl
#  Purpose:
#  In:
#  Out:
#


# --------------------------------------------------------------------------
# Sliding window analysis
# Written by Stephan Ossowski
# --------------------------------------------------------------------------


my $file = shift;
open FILE, $file or die;

my $win_pos = 50;
my $next_step = 60;
my $count = 0;
my $sum = 0;

while( <FILE> ) {
	chomp;
	my @a = split("\t", $_);

	if($a[0] > 50) {
		if($a[0] > $next_step) {
			print "$win_pos\t$count\t$sum\n";
			$win_pos += 10;
			$next_step += 10;
			$count = 0;
		}
		$count += $a[1];
		$sum += $a[1];
	}
}
close(FILE);
