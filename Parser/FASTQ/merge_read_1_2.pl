#!/usr/bin/perl

use strict;
use warnings;

my $usage  = "$0 file1 file2\n";
my $file1   = shift or die $usage;
my $file2   = shift or die $usage;

open F1, $file1 or die "Cannot open input read 1 file\n";
open F2, $file2 or die "Cannot open input read 2 file\n";

while( <F1> ) {
	my $sh1   = $_;
	my $seq1  = <F1>;
	my $qh1   = <F1>;
	my $qual1 = <F1>;

	my $sh2   = <F2>;
	my $seq2  = <F2>;
	my $qh2   = <F2>;
	my $qual2 = <F2>;

	print "$sh1$seq1$qh1$qual1$sh2$seq2$qh2$qual2";
}

close F1;
close F2;

exit(0);
