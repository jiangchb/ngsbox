#!/usr/bin/perl
######################################################################################
#Author 	Stephan Ossowski
#Date 		10/23/07
#Version	0.9
#Function	Filter left over reads using qvalue, chas and svalue
######################################################################################

use strict;
use warnings;

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
