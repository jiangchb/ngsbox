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
#  Module: Annotation::Genotyping::Genotyping_Add_Refbase.pl
#  Purpose:
#  In:
#  Out:
#


####################################################################
# Overlap Genotyping SNPs with NGS SNPs
# Author: Stephan Ossowski
####################################################################


my $usage = "\n\n$0 reference genotype_file\n\n";
my $reffile  = shift or die $usage;
my $genotype = shift or die $usage;

my %genotypes = ();


# Read reference genome
my %ref = ();
my $chr = "NA";
open REF, $reffile or die "Cannot open $reffile file\n";
while( <REF> ) {
	chomp;
	if( $_ =~ />/ ) {
		if($chr ne "NA") {
			print STDERR $chr . "\t" . substr($ref{$chr}, 0, 100) . "\n";
		}
		$chr = substr($_, 1);
	}
	else {
		$ref{$chr} .= $_;
	}
}
print STDERR $chr . "\t" . substr($ref{$chr}, 0, 100) . "\n";
close REF;


# Get all genotyped positions and add reference base
open GENO, $genotype or die "Cannot open $genotype file\n";
while( <GENO> ) {
	chomp;
	my ($sample, $id, $chr, $pos, $allele_1, $allele_2) = split("\t", $_);

	if($chr eq "X") { $chr = 23; }
	if($chr eq "Y") { $chr = 24; }
	if($chr eq "MT") { $chr = 25; }
	$allele_2 = substr($allele_2, 0, 1);

	if($chr ne "0" && $chr ne "XY" && $pos > 0) {
		if( length($ref{$chr}) > $pos-1 ) {
			my $refbase = substr( $ref{$chr}, $pos-1, 1);
			print "$sample\t$id\t$chr\t$pos\t$allele_1\t$allele_2\t$refbase\n";
		}
	}
}
close GENO;
exit(0);
