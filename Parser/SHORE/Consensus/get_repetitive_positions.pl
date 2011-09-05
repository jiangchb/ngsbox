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
#  Module: Parser::SHORE::Consensus::get_repetitive_positions.pl
#  Purpose:
#  In:
#  Out:
#


my $usage= "\n$0 consensus\n\n" ;
my $file = shift or die $usage;

open FILE, $file or die $usage;

while (<FILE>) {
	chomp;
	my @a = split("\t", $_);

	if( $a[10] > 1 ) {
		print $a[1] . "\t" . $a[10] . "\n";
	}
}
