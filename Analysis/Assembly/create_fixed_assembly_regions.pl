#!/usr/bin/perl

use strict;
use warnings;

my $block_size = shift;
my $offset     = shift;
my $max        = shift;

my $start = $offset;
my $end   = $offset + $block_size - 1;

while( $start <= $max ) {
	print "$start\t$end\n";
	$start += $block_size;
	$end += $block_size;
}

exit(0);
