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
#  Module: Analysis::ConsensusSequence::Callability::accessible_space_nonqual.pl
#  Purpose:
#  In:
#  Out:
#



my $chr = shift or die;
my $beg = shift or die;
my $end = shift or die;
my $sum = shift or die;
my $ref = shift or die;
my $snp = shift or die;
my $het = shift or die;
my $del = shift or die;


my $callable = 0;
my $uncovered = 0;
my $repetitive = 0;
my $other = 0;
my %good_hash = ();

print STDERR "Start REF parsing\n";

# Callable ref pos
open REF, $ref or die "Cannot open input file.\n";
while( <REF> ) {
	my @a = split("\t", $_);
	if( ($a[1] eq $chr) && ($a[2] >= $beg) && ($a[2] <= $end) ) {
		$callable++;
		$good_hash{$a[2]} = $a[1];
	}
}
close REF;

print STDERR "Start SNP parsing\n";

# Callable SNP pos
open SNP, $snp or die "Cannot open input file.\n";
while( <SNP> ) {
	my @a = split("\t", $_);
	if( ($a[1] eq $chr) && ($a[2] >= $beg) && ($a[2] <= $end) ) {
        	$callable++;
		$good_hash{$a[2]} = $a[1];
	}
}
close SNP;

print STDERR "Start HET parsing\n";

# Callable het pos
open HET, $het or die "Cannot open input file.\n";
while( <HET> ) {
	my @a = split("\t", $_);
	if( ($a[1] eq $chr) && ($a[2] >= $beg) && ($a[2] <= $end) ) {
		$callable++;
		$good_hash{$a[2]} = $a[1];
	}
}
close HET;

print STDERR "Start DEL parsing\n";

# Callable del pos
open DEL, $del or die "Cannot open input file.\n";
while( <DEL> ) {
	my @a = split("\t", $_);
	if( ($a[1] eq $chr) && ($a[2] >= $beg) && ($a[3] <= $end) ) {
		$callable++;
		$good_hash{$a[2]} = $a[1];
	}
}
close DEL;

print STDERR "Start SUM parsing\n";

# Other positions
open SUM, $sum or die "Cannot open input file.\n";
while( <SUM> ) {
	my @a = split("\t", $_);
	if( ($a[0] eq $chr) && ($a[1] >= $beg) && ($a[1] <= $end) ) {
		if( ! exists $good_hash{$a[1]} ) {
			if( ($a[3] - $a[9]) < 3 ) {
				$uncovered++;
			}
			elsif( $a[10] > 1 ) {
				$repetitive++;
			}
			else {
				$other++;
			}
		}
	}
}
close SUM;

print "$callable\t$uncovered\t$repetitive\t$other\n";

exit(0);
