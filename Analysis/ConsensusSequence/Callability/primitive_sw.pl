#!/usr/bin/perl

# --------------------------------------------------------------------------
# Sliding window analysis
# Written by Stephan Ossowski
# --------------------------------------------------------------------------

use strict;
use warnings;

my $file = shift;
open FILE, $file or die;

my $win_pos = 50;
my $next_step = 60;
my $count = 0;
my $sum = 0;

while( <FILE> ) {
	chomp;
	my @a = split("\t", $_);

	if($a[0] > 50) {
		if($a[0] > $next_step) {
			print "$win_pos\t$count\t$sum\n";
			$win_pos += 10;
			$next_step += 10;
			$count = 0;
		}
		$count += $a[1];
		$sum += $a[1];
	}
}
close(FILE);
