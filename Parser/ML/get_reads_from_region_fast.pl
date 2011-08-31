#! /usr/bin/perl
use warnings;
use strict;

my $map_list = shift;
my $chr      = shift;
my $start    = shift;
my $end      = shift;

### Get reads from target region
open FILE, $map_list or die "cannnot open $map_list\n";
while (<FILE>) {
	my @a = split " ", $_;

	if( $a[0] == $chr && $a[1] >= $start && $a[1] <= $end) {
		print $_;
	}

	if( $a[0] == $chr && $a[1] > $end) { exit(0); }
}


exit(0);
