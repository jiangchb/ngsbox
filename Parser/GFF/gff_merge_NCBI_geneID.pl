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
#  Module: Parser::GFF::
#  Purpose:
#  In:
#  Out:
#

###############################################################
#Author 	Stephan Ossowski, Korbinian Schneeberger 
#Date 		07/04/07
#Version	0.1
#Input		Read flat sequence file and write fasta 
#		formated file to stdout
###############################################################

my $usage = "\n$0 ncbi_file gff_file\n\n";
my $ncbi_file = shift or die $usage;
my $gff_file = shift or die $usage;

open NCBIFILE, $ncbi_file or die;
open GFFFILE, $gff_file or die;

my %rna = ();

while(<NCBIFILE>) {
	chomp;
	my @a = split("\t", $_);

	my $gene_id = $a[1];
	my $rna_id = $a[3];
	my ($gff_id, $junk) = split(/\./, $rna_id);

	$rna{$gff_id} = $gene_id;
}

my %trans = ();

while(<GFFFILE>) {
	chomp;
	my @a = split("\t", $_);
	my @b = split(";", $a[8]);
	my ($junk, $gff_id) = split(" ", $b[0]);
	$gff_id =~ s/"//g;

	if(exists $rna{$gff_id}) {
		$trans{$gff_id} = $gff_id ."\t". $rna{$gff_id} ."\t". $a[8];
	}
	else {
		$trans{$gff_id} = $gff_id ."\tNA\t". $a[8];
	}
}


foreach my $gff_id (keys %trans) {
	print $trans{$gff_id} . "\n";
}

exit(0);
