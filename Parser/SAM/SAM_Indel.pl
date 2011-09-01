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

my $file = shift or die "\n\nUsage: $0 SAMfile\n\n";
open FILE, $file;

while(<FILE>) {
	my @a = split("\t", $_);

	#my $pe = 0;
	#if ( hex($a[1]) & 0x0040 ) {
	#	$pe = 1;
	#}
	#elsif ( hex($a[1]) & 0x0080 ) {
	#	$pe = 2;
	#}
	
	if( ($a[5] =~ /D/) || ($a[5] =~ /I/) ) {
		my @c = split(/[MIDS]/, $a[5]);
		
		if( $c[0] >= 8 && $c[$#c] >= 8) {
			print $_;
		}
	}
}

close FILE;

exit(0);
