#!/usr/bin/perl

use strict;
use warnings;

my $file = shift or die;

my $snp = 0;
my $snp_ath = 0;
my $snp_lyr = 0;
my $no_snp = 0;

open FILE, $file or die "Cannot open input file.\n";
while( <FILE> ) {
	chomp;
	my @a = split("\t", $_);
	if($a[6] ne "NA" || $a[14] ne "NA") {
		$snp++;

		if($a[6] ne "NA") {
			$snp_ath++;
		}
		if($a[14] ne "NA") {
			$snp_lyr++;
		}
	}
	else {
		$no_snp++;
	}
}

my $total = $snp + $no_snp;

print "$total\t$no_snp\t$snp\t$snp_ath\t$snp_lyr\n";

exit(0);
