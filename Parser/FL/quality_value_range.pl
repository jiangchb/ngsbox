#! /usr/bin/perl

use strict;

my $usage = "perl filename column";
my $file = shift or die $usage;
my $column =shift or die $usage;
open FILE, $file;

my $count = 0;
my %Q = ();
while (my $line = <FILE> and $count < 10000) {
	my @a = split "\t", $line;
	for (my $i = 0; $i<length($a[$column-1]); $i++) {
		$Q{ord(substr($a[$column-1], $i, 1))}++;
	}
	$count++;
}


foreach my $key (sort {$a <=> $b} keys %Q) {
	print $key, "\t", $Q{$key}, "\n";
}

