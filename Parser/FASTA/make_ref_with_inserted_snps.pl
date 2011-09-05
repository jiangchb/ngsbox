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
#  Module: Parser::FASTA::make_ref_with_inserted_snps.pl
#  Purpose:
#  In:
#  Out:
#



# Makes a "reference sequence" with inserted snps.

my $usage = "make_ref_with_inserted_snps.pl snpfile fastafile\n";
my $snpfile = shift or die $usage;
my $fastafile = shift or die $usage;

open SNP, $snpfile;

my %SNPS = ();

while(my $line = <SNP>) {
	my @a = split " ", $line;
	$SNPS{$a[1]."#".$a[2]} = $a[4];
}

close SNPS;

open FASTA, $fastafile;

my $chr = 0;
my $pos = 0;
while (my $line = <FASTA>) {
	chomp($line);
	if (substr($line, 0, 1) eq ">") {
		$chr++;
		$pos = 0;
		print $line, "\n";
	}
	else {
		for (my $i = 0; $i < length($line); $i++) {
			$pos++;
			if (defined($SNPS{$chr."#".$pos})) {
				print $SNPS{$chr."#".$pos};
			}
			else {
				print substr($line, $i, 1);
			}
		}
		print "\n";
	}
}

