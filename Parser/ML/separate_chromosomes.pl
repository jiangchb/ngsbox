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
#  Module: Parser::ML::separate_chromosomes.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 maplist outfolder\n";

my $file      = shift or die $usage;
my $outfolder = shift or die $usage;

open FILE, $file or die "cannnot open $file\n";

my $chromosome = -1;

while(<FILE>) {
	my $line = $_;
	my @entries = split(/\t/, $line);
	if($entries[0] != $chromosome) {
		if($chromosome != -1) { close OUT; }
		$chromosome = $entries[0];
		open(OUT, ">$outfolder/map_chr$chromosome.list") or die "cannot open $outfolder/map_chr$chromosome.list\n";
	}

	print OUT $line;
}

exit(0);
