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
#  Module: Analysis::Assembly::Overlap::new_gene_finder_generic.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 gff segments\n";

my $file1 = shift or die $usage;
my $file2 = shift or die $usage;


### Read gene annotation file (gff)
my %anno = ();
my %anno_strand_spec = ();
open FILE1, $file1 or die "Cannot open annotation file\n";

while (<FILE1>) {
	my @a = split "\t";
	for (my $i = $a[3]; $i<= $a[4]; $i++) {

		$anno{$a[0]."#".$i} = 1;

		$anno_strand_spec{$a[0]."#".$a[6]."#".$i} = 1;
	}
}
close FILE1;


### Read coverage segment file and check overlap with annotation
open FILE2, $file2 or die "Cannot open segments file\n";

while ( my $line = <FILE2>) {

	chomp($line);
	my @a = split("\t", $line);
	my $overlap_anno = 0;
	my $overlap_strand = 0;
	
	for (my $i = $a[1]; $i <= ($a[1] + $a[2]); $i++) {
		
		if( exists $anno{$a[0]."#".$i} ) {
			$overlap_anno++;
		}
		if ( exists $anno_strand_spec{$a[0]."#".$a[3]."#".$i} ) {
			$overlap_strand++;
		}
	}
	
	if( $overlap_anno == 0 && $a[7] >= 7) {
		print "NewGene\t$line\n";
	}
	elsif( $overlap_strand == 0 && $a[7] >= 30 ) {
		print "NewAnti\t$line\n";
	}
}

exit;
