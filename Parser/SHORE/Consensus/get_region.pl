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
#  Module: Parser::SHORE::Consensus::get_region.pl
#  Purpose:
#  In:
#  Out:
#


my $usage= "\n$0 consensus_summary.txt chr begin end\n\n" ;
my $file  = shift or die $usage;
my $chr   = shift or die $usage;
my $begin = shift or die $usage;
my $end   = shift or die $usage;

open FILE, $file or die $usage;
my $flag = 0;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[0] == $chr and $a[1] >= $begin and $a[1] <= $end) {
		print $line;
		$flag = 1;
	}
	else {
		if ($flag == 1) {
			exit(1);
		}
	}
}

