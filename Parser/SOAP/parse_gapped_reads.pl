#! /usr/bin/perl

use strict;

my $usage = "$0 soapoutput\n";
my $file = shift or die $usage;
open FILE, $file or die $usage;


my %IDS = ();
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[9] > 99) {
		if (not defined($IDS{$a[0]})) {
			print $a[0], "\n";
			$IDS{$a[0]} = 1;
		}
	}
}


