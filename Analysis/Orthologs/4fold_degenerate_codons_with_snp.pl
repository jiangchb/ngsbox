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
#  Module: Analysis::Orthologs::4fold_degenerate_codons_with_snp.pl
#  Purpose:
#  In:
#  Out:
#



my $file = shift or die;

my $snp = 0;
my $snp_ath = 0;
my $snp_lyr = 0;
my $no_snp = 0;

open FILE, $file or die "Cannot open input file.\n";
while( <FILE> ) {
	chomp;
	my @a = split("\t", $_);
	if($a[6] ne "NA" || $a[14] ne "NA") {
		$snp++;

		if($a[6] ne "NA") {
			$snp_ath++;
		}
		if($a[14] ne "NA") {
			$snp_lyr++;
		}
	}
	else {
		$no_snp++;
	}
}

my $total = $snp + $no_snp;

print "$total\t$no_snp\t$snp\t$snp_ath\t$snp_lyr\n";

exit(0);
