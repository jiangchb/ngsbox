#!/usr/bin/perl

# select columns (comma seperated) and one constraint on one column (comma seperated "column,constrains") to filter files

use strict;
use warnings;

my $usage = "$0 selected_columns constraint file\n";
my $columns    = shift or die $usage;
my $constraint = shift or die $usage;
my $file       = shift or die $usage;

my @colum = split(",", $columns);
my ($ccolumn, $cvalue) = split(",", $constraint);

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;

	my @e = split("\t", $_);

	if( $e[$ccolumn - 1] > $cvalue ) {

		my $out = "";

		for(my $i = 0; $i <= $#colum; $i++) {

			my $value = $e[$colum[$i]-1];
		
			$out .= $value . "\t";
		}

		chop($out);

		print $out . "\n";
	}
}

close IN;

exit(0);
