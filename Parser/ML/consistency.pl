#!/usr/bin/perl
use strict;
use warnings;

my $file    = shift;
my $column_count = shift;
my $integer_columns = shift;
my $print = shift;

if( ! defined $print ) { $print = 0; }

if( ! defined $column_count ) { 
	print "\nUsage: consistency.pl <file> <column_count> <integer_columns starting with 0>\n\n" .
		"Example(maplist): consistency.pl map.list 13 0,1,3,5,6,7\n\n" .
		"Example (reads): consistency.pl reads_0.fl 6 0,2\n\n";
	exit(0);
}

my @ints = split( ",", $integer_columns );

my $counter = 0;
open IN, $file or die "Cannot open input file\n";

while ( <IN> ) {
	chomp;
	my $row = $_;
	my @elements = split(/\t/, $row);
	$counter++;

	my $correct = 1;
	if(@elements != $column_count) {
		print STDERR "Line:$counter, Wrong number of columns:\n$row\n";
		$correct = 0;
	}
	else {
		foreach my $column (@ints) {
			if( (! defined $elements[$column]) || ($elements[$column] !~ /\d+/) ) {
				print STDERR "Row:$counter, Column $column, Not a number:\n$row";
				$correct = 0;
			}
		}
	}

	if( ($correct == 1) && ($print == 1) ) {
		print "$row\n";
	}
}

close IN;

exit(0);
