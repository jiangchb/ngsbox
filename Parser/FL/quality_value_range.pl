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
#  Module: Parser::FL::quality_value_range.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "perl filename column";
my $file = shift or die $usage;
my $column =shift or die $usage;
open FILE, $file;

my $count = 0;
my %Q = ();
while (my $line = <FILE> and $count < 10000) {
	my @a = split "\t", $line;
	for (my $i = 0; $i<length($a[$column-1]); $i++) {
		$Q{ord(substr($a[$column-1], $i, 1))}++;
	}
	$count++;
}


foreach my $key (sort {$a <=> $b} keys %Q) {
	print $key, "\t", $Q{$key}, "\n";
}

