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
#  Module: Analysis::Tiling_array::correct_probes_4_polymorphism.pl
#  Purpose:
#  In:
#  Out:
#



my $MIN_BASE_QUAL = 25;

my $usage = "$0 tilingprobes max_core_missing max_boarder_missing quality_reference_file [...]\n";

my $file = shift or die $usage;
my $MAX_MM_CORE = shift;
my $MAX_MM_BOARDER = shift;

my $ACC_NUM = @ARGV + 0;

my %CONS = ();

for (my $i = 0; $i <= $#ARGV; $i++) { 
	open REF, $ARGV[$i] or die $usage;
	print STDERR "Open ", $ARGV[$i], "\n";
	while (<REF>) {
		my @a = split " ";
		if ($a[5] >= $MIN_BASE_QUAL) {
			$CONS{$a[1]."#".$a[2]}++;
		}
		if ($a[2] % 1000000 == 0) {
			print STDERR $a[1], "\t", $a[2], "\n";
		}
	}
	close REF;
}
 
print STDERR "Finished parsing conserved bases\n";

open FILE, $file or die $usage;
open OUT, "> ".$file.".core$MAX_MM_CORE.boarder$MAX_MM_BOARDER.txt";

while (my $line = <FILE>) {
	if (substr($line, 0, 1) ne "#") {
		my @a = split " ", $line;
		my $mm_border = 0;
		my $mm_core= 0;
		for (my $i = $a[4]; $i < $a[4]+2; $i++) {
			if (not defined($CONS{$a[3]."#".$i}) or $CONS{$a[3]."#".$i} < $ACC_NUM) {
				$mm_border++;
			}
		}
		for (my $i = $a[4]+2; $i <= $a[5]-2; $i++) {
			if (not defined($CONS{$a[3]."#".$i}) or $CONS{$a[3]."#".$i} < $ACC_NUM) {
                                $mm_core++;
                        }
		}
		for (my $i = $a[5]-1; $i<=$a[5]; $i++) {
			if (not defined($CONS{$a[3]."#".$i}) or $CONS{$a[3]."#".$i} < $ACC_NUM) {
                                $mm_border++;
                        }
		}

		print STDERR $a[3], "\t", $a[4], "\t", $mm_core, "\t", $mm_border, "\n";

		if ($mm_core <= $MAX_MM_CORE) {
			if ($mm_border <= $MAX_MM_BOARDER) {
				print OUT $line;
			}
		}
		
	}
}

