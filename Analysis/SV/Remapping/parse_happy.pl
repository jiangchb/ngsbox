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
#  Module: Analysis::SV::Remapping::parse_happy.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 maplist\n";
my $file = shift or die $usage;

open FILE, $file or die $usage;

my %IDS = ();

my $count = 0;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	$count++;
	print STDERR $a[0], "\t", $a[1], "\n" if $count%1000000 == 0;
	if ($a[9] != 0 and $a[9]%3 == 0) {
		$IDS{$a[3]} = 1;
	}
}

foreach my $id (keys %IDS) {
	print $id, "\n";
}




