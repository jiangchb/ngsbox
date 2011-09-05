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
#  Module: Analysis::Validation::PCR::pcr_primer3_bacteria_refcall.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $usage = "\n$0 type ref_calls snp_calls\n\n";
my $type = shift or die $usage;
my $ref  = shift or die $usage;
my $snp  = shift or die $usage;

my %refcall = ();


### Read reference calls
open REF, $ref or die "Cannot open $ref\n";
while(<REF>) {
	my @a = split(/\t/, $_);
	$refcall{$a[1]}{$a[2]} = $a[3] . ";" . $a[8];
}


### Read SNPs
open SNP, $snp or die "Cannot open $snp\n";
my $target_start = 451;
my $target_length = 100;

while(<SNP>) {
        chomp;
	my @a = split(/\t/, $_);
	my $chr = $a[1];
	my $validation_start = $a[2] - 500;
	my $Xed_beg = $a[2];
	my $Xed_end = $a[2];
	if($type eq "indel") {
		$Xed_end = $a[3];
	}

	# Concatenate reference calls around focal SNP, rest 'N'
	my $seq = "";
	my $ambi = 0;

	for(my $i = $validation_start; $i < $validation_start + 1001; $i++) {
		
		if( $i >= $Xed_beg && $i <= $Xed_end ) {
			$seq .= "X";
		}
		elsif( exists $refcall{$chr}{$i} ) {
			my ($base, $avg_hits) = split(/;/ , $refcall{$chr}{$i});
			
			if($avg_hits < 1.2) {
				$seq .= $base;
			}
			else {
				$seq .= "N";
				$ambi++;
			}
		}
		else {
			$seq .= "N";
			$ambi++;
		}
	}
	#if( $ambi < 100 ) {
		if( $type eq "snp") {
			print ">$chr-$a[2]-$a[3]$a[4]($validation_start):TARGET=$target_start,$target_length\n";
		}
		elsif( $type eq "indel" ) {
			print ">$chr-$a[2]-Indel$a[5]($validation_start):TARGET=$target_start,$target_length\n";
		}
		print "$seq\n";
	#}
}
