#!/usr/bin/perl
###########################################################################
# Author 	Korbinian Schneeberger 
# Date 		09/26/07
# Version	0.2
# Input		map.list, map.list.idx, chr, start, end
# Function	returns all reads from a specified region of the mapping
###########################################################################

use strict;
use warnings;

use lib "/ebio/abt6/stephan/pgsp/Assembly";
use parse_mapping;

my $file	= shift;
my $chr		= shift;
my $begin	= shift;
my $end		= shift;
my $read_length	= shift;
my $out_dir	= shift;


for(my $i = $begin; $i< $end; $i+=100) {
	my $end = $i + 200;
	my @results = parse_mapping::get($file, $chr, $i - $read_length + 1, $end);

	open OUT, ">$out_dir/$chr-$i-$end.fa" or die;

	foreach (@results) {
		my @entries = split(" ", $_);
		my $seq = $entries[2];
	
		print OUT ">$entries[3] | $entries[0] | $entries[1]\n";

		for(my $i = 0; $i < length($seq); $i++) {
			if(substr($seq, $i, 1) eq "[") {
				$i+=2;
				if(substr($seq, $i, 1) ne "-") {
					print OUT substr($seq, $i, 1);
				}
				$i++;
			}
			else {
				if(substr($seq, $i, 1) ne "-") {
					print OUT substr($seq, $i, 1);
				}
			}
		}
		print OUT "\n";
	}

	close OUT;
}

exit(0);

