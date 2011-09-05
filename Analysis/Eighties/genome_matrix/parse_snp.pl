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
#  Module: Analysis::Eighties::genome_matrix::parse_snp.pl
#  Purpose:
#  In:
#  Out:
#

###########################################################################
# parse positions having at least one SNP from genome_matrix.base_call file 
###########################################################################



my $usage = "$0 genomeMatrix_basse_call file\n";

my $file = shift or die $usage;

open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;


	for (my $i = 3; $i < @a; $i+=1) {
		if ($a[$i] ne $a[2] and $a[$i] ne "N") {
	

	chomp($line);
	print $line, "\n";
        last;
        } 
          
        }
}


