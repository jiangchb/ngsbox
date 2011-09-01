#!/usr/bin/perl
###############################################################
#Author 	Stephan Ossowski, Korbinian Schneeberger 
#Date 		07/04/07
#Version	0.1
#Input		Read flat sequence file and write fasta 
#		formated file to stdout
###############################################################

use strict;
use warnings;

my $file = shift;
open FILE, $file;

while(<FILE>) {
	my @a = split " ";

	if($a[2] eq "CDS") {
		my $len = $a[4] - $a[3] + 1;
		print "1\t" . $a[3] . "\t" . $a[4] . "\t" . $a[6] . "\n";
	}
}

exit(0);
