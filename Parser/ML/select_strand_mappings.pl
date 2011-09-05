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
#  Module: Parser::ML::select_strand_mappings.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "perl $0 maplist D|P";

my $file = shift or die $usage;
my $dir = shift or die $usage;
die $usage if ($dir ne "P" and $dir ne "D");

open FILE, $file or die $usage;

while (<FILE>) {
	my @a = split " ";
	if ($a[4] eq $dir) {
		print $_;
	}
}


