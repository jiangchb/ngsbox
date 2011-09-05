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
#  Module: Analysis::Validation::PCR::pcr_primer3_bacteria_reffasta.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $usage = "\n$0 snp_calls ref_fasta\n\n";
my $snp   = shift or die $usage;
my $fasta = shift or die $usage;

my $refseq  = ();
my %refcall = ();


### Read fasta sequence
open FASTA, $fasta or die "Cannot open $fasta\n";
while(<FASTA>) {
	chomp;
	if($_ !~ />/) { 
		$refseq .= $_; 
	}
}

open SNP, $snp or die "Cannot open $snp\n";
my $target_start = 450;
my $target_length = 100;

while(<SNP>) {
        chomp;
	my @a = split(/\t/, $_);
	my $validation_start = $a[2] - 500;

	# Use original reference sequence for primer design
	print ">$a[0]-$a[1]-$a[2]-$a[3]($a[4]):TARGET=$target_start,$target_length\n";
	print substr($refseq, $validation_start, 800)  . "\n";
}
