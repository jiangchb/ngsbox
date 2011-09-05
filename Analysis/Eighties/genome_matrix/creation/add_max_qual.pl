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
#  Module: Analysis::Eighties::genome_matrix::creation::add_max_qual.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 genomeMatrix\n";

my $file = shift or die $usage;

open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;

	my $max_a = -1;
	my $max_c = -1;
	my $max_g = -1;
	my $max_t = -1;
	my $max_d = -1;

	my $freq_a = 0;
	my $freq_c = 0;
	my $freq_g = 0;
	my $freq_t = 0;
	my $freq_d = 0;
	my $freq_n = 0;

	for (my $i = 3; $i < @a; $i+=2) {
		if ($a[$i] eq "A") {
			if ($a[$i+1] > $max_a) {
				$max_a = $a[$i+1];
			}
			$freq_a++;
		}
		elsif ($a[$i] eq "C") {
			if ($a[$i+1] > $max_c) {
                                $max_c = $a[$i+1];
                        }
                        $freq_c++;
		}
		elsif ($a[$i] eq "T") {
			if ($a[$i+1] > $max_t) {
                                $max_t = $a[$i+1];
                        }
                        $freq_t++;
                }
		elsif ($a[$i] eq "G") {
			if ($a[$i+1] > $max_g) {
                                $max_g = $a[$i+1];
                        }
                        $freq_g++;
                }
		elsif ($a[$i] eq "-") {
                        if ($a[$i+1] > $max_d) {
                                $max_d = $a[$i+1];
                        }
                        $freq_d++;
                }
		elsif ($a[$i] eq "N") {
                        $freq_n++;
                }
	}

	chomp($line);
	print $line, "\t", $freq_a, "\t", $max_a, "\t", $freq_c, "\t", $max_c, "\t", $freq_g, "\t", $max_g, "\t", $freq_t, "\t", $max_t, "\t", $freq_d, "\t", $max_d, "\n";

}


