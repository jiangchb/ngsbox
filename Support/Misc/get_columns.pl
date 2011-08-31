#!/usr/bin/perl

use strict;
use warnings;

my $usage = "$0 selected_columns file\n";
my $selected_columns = shift or die $usage;
my $file             = shift or die $usage;


my @colums = split(",", $selected_columns);


open IN, $file or die "Cannot open input file\n";

my $counter = 1000000;

while( <IN> ) {
	chomp;

	my @e = split("\t", $_);

	my $out = "";

	for(my $i = 0; $i <= $#colums; $i++) {

		my $current = $e[$colums[$i]-1];

		# Optimization for human chromosomes
		if(substr($current, 0, 3) eq "chr") {
			$current = substr($current, 3);
			if($current eq "X") {
				$current = 23;
			}
			if($current eq "Y") {
				$current = 24;
			}
			if($current eq "M") {
				$current = 25;
			}
		}

		$out .= $current . "\t";
	}

	chop($out);

	print $counter . "_" . $out . "\n";
	$counter++;
}

close IN;

exit(0);
