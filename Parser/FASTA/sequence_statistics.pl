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
#  Module: Parser::FASTA::calc_assem_stats.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 threshold fasta\n";
my $threshold = shift;
my $file = shift or die $usage;

###########################################################
# Read in file
open FILE, $file or die $usage;

my %CTG_LTH = ();
my $seq = "";

while (my $line = <FILE>) {
	chomp($line);
	if (substr($line, 0, 1) eq ">") {
		if ($seq ne "") {
			$CTG_LTH{length($seq)}++;
		}
		$seq = "";
	}
	else {
		$seq .= $line;
	}
}
if ($seq ne "") {
	$CTG_LTH{length($seq)}++;
}


###########################################################
# Analyze lengths and number
my $max = 0;
my $min = -1;
my $total_count = 0;
my $count_above_threshold = 0;
my $good_nucleotides = 0;

foreach my $length (sort {$b <=> $a} keys %CTG_LTH) {
	
	# set max length
	if ($max < $length) {
		$max = $length;
	}

	# set min length
	if ($min == -1 || $min > $length) {
		$min = $length;
	}

	$total_count += $CTG_LTH{$length};
	if ($length >= $threshold) {
		$count_above_threshold += $CTG_LTH{$length};
		$good_nucleotides += ($length * $CTG_LTH{$length});
	}
}

print "Reads:\t\t$total_count\n";
print "Reads above threshold length:\t\t$count_above_threshold\n";
print "Good nucelotides:\t\t$good_nucleotides\n";
print "Longest read:\t\t$max\n";
print "Shortest read:\t\t$min\n";


