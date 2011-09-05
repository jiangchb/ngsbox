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
#  Module: Analysis::Eighties::genome_matrix::creation::genome_matrix_frequency_calling.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 genomeMatrixMaxQual MinMaxQual MinQual\n";

my $file = shift or die $usage;
my $MinMaxQual = shift or die $usage;
my $MinQual = shift or die $usage;

open FILE, $file or die $usage;
open OUT, "> genome_matrix.frequency.$MinMaxQual.$MinQual.txt";

while (my $line = <FILE>) {
	my @a = split " ", $line;

	my $max_a = $a[164];
	my $max_c = $a[166];
	my $max_g = $a[168];
	my $max_t = $a[170];
	my $max_d = $a[172];

	print OUT $a[0], "\t", $a[1], "\t", $a[2];

	my $count_a_high = 0;
	my $count_a_low = 0;
	my $count_c_high = 0;
        my $count_c_low = 0;
	my $count_g_high = 0;
        my $count_g_low = 0;
	my $count_t_high = 0;
        my $count_t_low = 0;
	my $count_d_high = 0;
        my $count_d_low = 0;
        my $count_n = 0;

	for (my $i = 3; $i <= 162; $i+=2) {
		if ($a[$i] eq "A") {
			if (($max_a >= $MinMaxQual and $a[$i+1] > $MinQual) or ($a[$i+1] > $MinMaxQual)) {
				$count_a_high++;
			}
			$count_a_low++;
		}
		elsif ($a[$i] eq "C") {
			if (($max_c >= $MinMaxQual and $a[$i+1] > $MinQual) or ($a[$i+1] > $MinMaxQual)) {
	                        $count_c_high++;
			}
                        $count_c_low++;
                }
		elsif ($a[$i] eq "G") {
			if (($max_g >= $MinMaxQual and $a[$i+1] > $MinQual) or ($a[$i+1] > $MinMaxQual)) {
	                        $count_g_high++;
	                }
               	        $count_g_low++;
                }
		elsif ($a[$i] eq "T") {
			if (($max_t >= $MinMaxQual and $a[$i+1] > $MinQual) or ($a[$i+1] > $MinMaxQual)) {
	                        $count_t_high++;
	                }
               	        $count_t_low++;
                }
		elsif ($a[$i] eq "-") {
                        if (($max_d >= $MinMaxQual and $a[$i+1] > $MinQual) or ($a[$i+1] > $MinMaxQual)) {
                                $count_d_high++;
                        }
                        $count_d_low++;
                }
		else {
			$count_n++;;
		}
	}
	
	print OUT "\t$count_a_high\t$count_c_high\t$count_g_high\t$count_t_high\t$count_d_high";
	print OUT "\t$count_a_low\t$count_c_low\t$count_g_low\t$count_t_low\t$count_d_low";
	print OUT "\t$count_n\n";
		

}


