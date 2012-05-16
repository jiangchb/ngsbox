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
#  Module: Parser::SHORE::Consensus::get_subset.pl
#  Purpose:
#  In:
#  Out:
#


my $usage= "\n$0 markerfile < consensus_summary.txt\n\n" ;
my $subset = shift or die $usage;

my %POS = ();

open FILE, $subset or die $usage;
while (my $line = <FILE>) {
        my @a = split " ", $line;
	my $id = $a[1]."#".$a[2];
	$POS{$id} = 1;
}
close FILE;

while (my $line = <STDIN>) {
        my @a = split " ", $line;
	my $id = $a[0]."#".$a[1];
	if (defined($POS{$id})) {
		print $line;
	}
}

