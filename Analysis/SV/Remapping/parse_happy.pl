#! /usr/bin/perl

use strict;

my $usage = "$0 maplist\n";
my $file = shift or die $usage;

open FILE, $file or die $usage;

my %IDS = ();

my $count = 0;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	$count++;
	print STDERR $a[0], "\t", $a[1], "\n" if $count%1000000 == 0;
	if ($a[9] != 0 and $a[9]%3 == 0) {
		$IDS{$a[3]} = 1;
	}
}

foreach my $id (keys %IDS) {
	print $id, "\n";
}




