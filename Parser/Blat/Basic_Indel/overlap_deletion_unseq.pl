#! /usr/bin/perl
use strict;

my $usage = "$0 deletions unsequenced\n";


my $del = shift or die $usage;
open UNSEQ, shift or die $usage;

my %U = ();

while (<UNSEQ>) {
	my @a = split " ";
	for (my $i = $a[2]; $i < $a[3]; $i++) {
		$U{$a[1]."#".$i} = 1;
	}
}

close UNSEQ;

open FILE, $del or die "Cannot open file\n";

while (<FILE>) {
	chomp;
	my @a = split " ";
	my $count = 0;
	
	for (my $i = $a[2]; $i <= $a[3]; $i++) {
		if (defined($U{$a[1]."#".$i})) {
	                $count++;
		}
        }

	print $_, "\t", $count, "\n";
}

