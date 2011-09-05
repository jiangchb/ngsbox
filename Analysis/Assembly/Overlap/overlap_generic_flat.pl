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
#  Module: Analysis::Assembly::Overlap::overlap_generic_flat.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 file1 file2\n";

my $file1 = shift or die $usage;
my $file2 = shift or die $usage;


### Read first result file
my %F1=();
open FILE1, $file1 or die "Cannot open file 1\n";

while (<FILE1>) {
	my @a = split " ";
	for (my $i = $a[2]; $i<= $a[3]; $i++) {
		$F1{$a[1]."#".$i} = 1;
	}
}
close FILE1;


### Read second result file and calculate overlap
my $total_count = 0;
my $total_length = 0;
open FILE2, $file2 or die "Cannot open file 2\n";

while (<FILE2>) {

	chomp();

	my @a = split " ";
	my $count = 0;
	my $length = 0;
	
	if( ($a[3] - $a[2] + 1) > 4 ) {
		for (my $i = $a[2]; $i<= $a[3]; $i++) {

			$length++;
			$total_length++;

			if (defined($F1{$a[1]."#".$i})) {
				$count++;
				$total_count++;
			}
		}

		my $overlap = $count / $length;

		print "$a[0]\t$a[1]\t$count\t$length\t$overlap\n";
	}
}

my $total_overlap = $total_count / $total_length;
print "$total_count\t$total_length\t$total_overlap\n";

exit;
