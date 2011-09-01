#!/usr/bin/perl

use strict;
use warnings;

my $split_num = shift;
my $file = shift;

my $count = 0;
my $file_count = 1;

open IN, $file or die "Cannot open input file\n";
open OUT, ">$file_count.txt" or die "Cannot open output file\n";

while( <IN> ) {
	print OUT "$_";
	$count++;

	if($count > $split_num) {
		$count = 0;
		$file_count++;
		close OUT;
		open OUT, ">$file_count.txt" or die "Cannot open output file\n";
	}
}
close IN; close OUT;

exit(0);
