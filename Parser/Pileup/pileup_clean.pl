#!/usr/bin/perl
###############################################################
#Author 	Stephan Ossowski, Korbinian Schneeberger 
#Date 		07/04/07
#Version	0.1
#Input		Check pileup variant file for correctness
###############################################################

use strict;
use warnings;

my $Usage = "\n\n$0 quality minlen minsup minfreq pileup_file\n\n";

my $quality = shift or die $Usage;
my $minlen  = shift or die $Usage;
my $minsup  = shift or die $Usage;
my $minfreq = shift or die $Usage;
my $file    = shift or die $Usage;

open FILE, $file;

while(<FILE>) {
	my @a = split("\t", $_);

	if( $#a >= 10 ) {
		if( ($a[5] >= $quality) && (length($a[8]) >= $minlen) && ($a[10] >= $minsup) && ($a[10]/$a[7] >= $minfreq) ) {
			print $_;
		}
	}
}

close FILE;

exit(0);
