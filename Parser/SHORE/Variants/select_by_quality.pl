#! /usr/bin/perl
use strict;
use warnings;

my $usage = "\n$0 quality_variants.txt qthreshold repetitive_threshold support frequency_lb frequency_ub snpsonly\n\n";
my $file = shift or die $usage;
my $qthres = shift or die $usage;
my $rthres = shift or die $usage;
my $support = shift;
my $frequency_lb = shift;
my $frequency_ub = shift;
my $snp = shift;

open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (	($a[8] <= $rthres) and 
		($a[5] >= $qthres) and 
		($snp != 1 || $a[4] ne "-") and 
		($a[7] >= $frequency_lb) and
		($a[7] <= $frequency_ub) and
		($a[6] >= $support)
	) {
		print $line;
	}
}



