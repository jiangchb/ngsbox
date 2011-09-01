#! /usr/bin/perl

use strict;

my $usage = "$0 gmoutput\n";
my $file = shift or die $usage;
open FILE, $file or die $usage;


my %IDS = ();
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[2] =~ m/\-/) {
		if (not defined($IDS{$a[3]})) {
			print $a[3], "\n";
			$IDS{$a[3]} = 1;
		}
	}
}


