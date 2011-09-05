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
#  Module: Analysis::Filter::read_quality_filter_fastq.pl
#  Purpose:
#  In:
#  Out:
#

######################################################################################
#Author 	Stephan Ossowski
#Date 		10/23/07
#Version	0.9
#Function	Filter left over reads using qvalue, chas and svalue
######################################################################################


my $usage = "\n\n$0 max_low_quality qvalue_threshold fastq-file\n\n";

my $max_low_quality  = shift or die $usage;
my $qvalue_threshold = shift or die $usage;
my $fq_file          = shift or die $usage;

open (FQFILE, $fq_file) or die "cannot open $fq_file\n";

### Loop through original reads and get qvalue entry
while( <FQFILE> ) {
	my $id1 = $_;
	my $seq1 = <FQFILE>;
	my $spacer1 = <FQFILE>;
	my $q1 = <FQFILE>;
	chomp($q1);

	my $id2 = <FQFILE>;
	my $seq2 = <FQFILE>;
	my $spacer2 = <FQFILE>;
	my $q2 = <FQFILE>;
	chomp($q2);


	my @qvalues1 = split(//, $q1);
	my @qvalues2 = split(//, $q2);


	for(my $i = 0; $i < @qvalues1; $i++) {
		$qvalues1[$i]   = ord($qvalues1[$i]) - 33;
	}
	for(my $i = 0; $i < @qvalues2; $i++) {
		$qvalues2[$i]   = ord($qvalues2[$i]) - 33;
	}

	my $low_quality_bases1  = 0;
	my $low_quality_bases2  = 0;

	for(my $i = 0; $i < @qvalues1; $i++) {
		if($qvalues1[$i] < $qvalue_threshold) {
			$low_quality_bases1++;
		}
	}

	for(my $i = 0; $i < @qvalues2; $i++) {
		if($qvalues2[$i] < $qvalue_threshold) {
			$low_quality_bases2++;
		}
	}

	
	if( ( $low_quality_bases1 <= $max_low_quality ) && ( $low_quality_bases2 <= $max_low_quality ) ) {
		print "$id1$seq1$spacer1$q1\n$id2$seq2$spacer2$q2\n";
	}
}

exit(0)
