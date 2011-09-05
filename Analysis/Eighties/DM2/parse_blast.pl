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
#  Module: Analysis::Eighties::DM2::parse_blast.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 blnfile\n";

my $file = shift or die $usage;
open FILE, $file or die $usage;

my $query = "";
my $max = -1;
my $min = -1;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (substr($line, 0, 6) eq "Query=") {
		$query = $a[1];
	}
	elsif (substr($line, 0, 8) eq " Score =" || substr($line, 0, 14) eq "Gap Penalties:") {
		if ($min != -1) {
			if ($max - $min > 300) {
				print $query, "\t", $min, "\t", $max, "\n";
			}
			$min = -1;
			$max = -1;
		}
	}
	elsif (substr($line, 0, 6) eq "Query:") {
		if ($min == -1 || $a[1] < $min) {
			$min = $a[1];
		}
		if ($max == -1 || $a[1] > $max) {
			$max = $a[1];
		}
		if ($min == -1 || $a[3] < $min) {
                        $min = $a[3];
                }
                if ($max == -1 || $a[3] > $max) {
                        $max = $a[3];
                }
	}
}


