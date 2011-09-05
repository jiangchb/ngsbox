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
#  Module: Parser::ML::GraphMapping::graphmap_sort4branchpos.pl
#  Purpose:
#  In:
#  Out:
#



# Sorting the map.list files at the 23.12 positions correctly:
# This:
# 
# 20.1
# 20.10
# 20.2
# 
# will be:
# 
# 20.1
# 20.2
# 20.10
#



my $usage = "$0 map.list\n";
my $file = shift or die $usage;

open FILE, $file or die "Cannot open file\n";

my $flag = 0;
my %LINE = ();

while (my $line = <FILE>) {
	my ($sample, $chr, $pos) = split " ", $line;
	if ($pos =~ /\./) {
		my ($ref, $strain) = split '\.', $pos;
		$LINE{$strain} .= $line;
		$flag = 1;
	}
	else {
		if ($flag == 1) {
			foreach my $strain (sort{$a <=> $b} keys %LINE) {
				print $LINE{$strain};
			}
			%LINE = ();
		}
		print $line;
		$flag = 0;
	}
}
if ($flag == 1) {
	foreach my $strain (sort{$a <=> $b} keys %LINE) {
        	print $LINE{$strain};
	}
        %LINE = ();
}



