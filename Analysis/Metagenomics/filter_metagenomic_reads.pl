#!/usr/bin/perl
######################################################################################
#Author 	Stephan Ossowski
#Date 		10/23/07
#Version	0.9
#Function	Filter left over reads using qvalue, chas and svalue
######################################################################################

use strict;
use warnings;

my $max_low_quality = shift;
my $qvalue_threshold = shift;
my $fl_file = shift;

open (FLFILE, $fl_file) or die "cannot open $fl_file\n";
open (ALL, ">$fl_file.filtered") or die "cannot open file\n";

### Loop through original reads and get qvalue entry
while( <FLFILE> ) {
	chomp($_);
	my ($id, $seq, $pe_flag, $qvalue, $svalue, $chas) = split(/\t/, $_);
	my @qvalues = split(//, $qvalue);
	my @chasts = split(//, $chas);

	for(my $i = 0; $i < @qvalues; $i++) {
		$qvalues[$i]   = ord($qvalues[$i]) - 64;
		$chasts[$i] = ord($chasts[$i]) + 10;
	}

	my $low_quality_bases  = 0;
	my $solexa_chas    = 0;


	for(my $i = 0; $i < @qvalues; $i++) {
		if( ($qvalues[$i] < $qvalue_threshold) || ($chasts[$i] < 57) ) {
			$low_quality_bases++;
		}
	}
	for(my $i = 0; $i < 12; $i++) {
		if($chasts[$i] < 60) { $solexa_chas++; }
	}

	
	if( ($low_quality_bases <= $max_low_quality) && ($solexa_chas < 2) && (&seq_complexity($seq) == 1) ) {
		print ALL "$id\t$seq\t$pe_flag\t$qvalue\t$svalue\t$chas\n";
	}
}

exit(0);


### Check complexity of reads
sub seq_complexity
{
        my $seq = shift;
        my $length = length($seq);
        my %occ   = (A => 0, C => 0,G => 0,T => 0);
        my %count = (A => 0, C => 0,G => 0,T => 0);

        for(my $i = 0; $i < length($seq); $i++) {
                $count{substr($seq, $i, 1)}++;
                $occ{substr($seq, $i, 1)} = 1;
        }

        if(
                ( ($occ{A} + $occ{C} + $occ{G} + $occ{T}) >= 2) &&
                ( ($length - $count{A}) > 6 ) &&
                ( ($length - $count{C}) > 6 ) &&
                ( ($length - $count{G}) > 6 ) &&
                ( ($length - $count{T}) > 6 )
        ) {
                return(1);
        }
        else {
                return(0);
        }
}
