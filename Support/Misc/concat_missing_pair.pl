#!/usr/bin/perl
use strict;
use warnings;

my $win_size = shift;
my $conc_size = shift;
my $file = shift;
open FILE, $file or die "Cannot open infile\n";

my $current_chr = -1;
my $start       = -1;
my $end         = -1;

while(<FILE>) {
	my ($chr, $pos, $count, $freq) = split " ";

	if( ($current_chr != $chr) || ($pos > $start + $conc_size) ) {
		my $length = $end - $start + 1;
		print "$current_chr\t$start\t$end\t$length\n";
		$current_chr = $chr;
		$start = $pos;
	}
	$end = $pos + $win_size;
}

exit(0);
