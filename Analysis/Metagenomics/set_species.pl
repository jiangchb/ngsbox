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
#  Module: Analysis::Metagenomics::set_species.pl
#  Purpose:
#  In:
#  Out:
#


my $mapfile     = shift;
my $speciesfile = shift;


### Read species names
my $counter = 1;
my %species = ();
open IN, $speciesfile or die "Cannot open species file\n";
while ( <IN> ) {
	chomp;
	$species{$counter} = $_;
	$counter++;
}
close IN;


### Set species name in mapfile
open IN, $mapfile or die "Cannot open map file\n";
while ( <IN> ) {
	chomp;
	my @elem = split(/\t/, $_);

	print $species{$elem[0]} . "\t$elem[1]\t$elem[2]\t$elem[3]\t$elem[4]\t$elem[5]\t$elem[6]\n";
}
close IN;

exit(0);
