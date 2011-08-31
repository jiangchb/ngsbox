#!/usr/bin/perl

###############################################################
# Convert Shore read files into fasta format
# Authors 	Stephan Ossowski, Korbinian Schneeberger 
# Date 		07/04/07
# Version	?
# Input		Shore sequence read file (fl) 
# Output	FASTA file to STDOUT
###############################################################

use strict;
use warnings;

my $file = shift;
open FILE, $file;

while(<FILE>) {
	my @a = split " ";
	print ">".$a[0]."\n".$a[1]."\n";

}

exit(0);
