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
#  Module: Parser::VCF::vcf_correction::fix_indel_vcf.pl
#  Purpose:
#  In:
#  Out:
#



my $vcf = shift or die "Please specify VCF file\n";

open VCF, $vcf or die "Cannot open input file\n";

# Loop through original reads and get prb entry
while( <VCF> ) {
	chomp;

	if($_ =~ /^#/) {
		print "$_\n";
	}
	else {
		my @a = split("\t", $_);
		my @b = split(";", $a[7]);

		my ($junk, $freq) = split("=", $b[2]);

		$a[7] =~ s/AC/DP/g;

		if($a[2] =~ /ins/) {
			$a[1]--;
		}

		if($freq > 0.8) {
			print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t" . "GT\t" . "1/1\n";
		}
		else {
			print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t" . "GT\t" . "0/1\n";
		}
	}
}

close VCF;

exit(0);
