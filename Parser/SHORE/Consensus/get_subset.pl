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


my $usage= "\n$0 consensus_summary.txt markerfile\n\n" ;
my $file = shift or die $usage;
my $subset = shift or die $usage;

#print STDERR "\n... will chop \"Chr\" of the chromosome identifier of consensus_summary\n";

my %POS = ();

open FILE, $subset or die $usage;
while (my $line = <FILE>) {
        my @a = split " ", $line;
	my $id = $a[1]."#".$a[2];
	$POS{$id} = 1;
	#print STDERR $id, "\n";
}
close FILE;

open FILE, $file or die $usage;
while (my $line = <FILE>) {
        my @a = split " ", $line;
	#if (substr($a[0], 0, 3) eq "Chr") {
	#	$a[0] = substr($a[0], 3);
	#}
	my $id = $a[0]."#".$a[1];
	if (defined($POS{$id})) {
		print $line;
	}
}

