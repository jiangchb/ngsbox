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
#  Module: Analysis::Validation::PCR::pcr_primer3_format.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $flat  = shift;
my $fasta = shift;

my $id = "";
my %contigs= ();

open FASTA, $fasta or die "Cannot open $fasta\n";
while(<FASTA>) {
	chomp;
	if($_ =~ />/) { $id = $_; }
	else { $contigs{$id} = $_; }
}

open FLAT, $flat or die "Cannot open $flat\n";
while(<FLAT>) {
        chomp;
	my @elem = split(/\t/, $_);
	my @snps = split(/;/, $elem[4]);
	my $target_start;
	my $target_length;

	if(scalar(@snps) == 1) {
		$target_start = $snps[0] - $elem[2] - 50;
		$target_length = 100;
	}
	elsif(scalar(@snps) > 1) {
		$target_start = $snps[0] - $elem[2] - 50;
		$target_length = $snps[-1] - $snps[0] + 100;
		if( 	($target_start < 30) ||
			($target_length > 400) || 
			($target_length < 0) 
		) { print STDERR "ERROR:$elem[0]-$elem[2]-$elem[3]\n"; }
	}
	else { print STDERR "ERROR\n"; }

	print ">$elem[0]-$elem[1]-$elem[2]-$elem[3]($elem[4]):TARGET=$target_start,$target_length\n";
	print $contigs{">$elem[1]-$elem[2]-$elem[3]"} . "\n";
}
