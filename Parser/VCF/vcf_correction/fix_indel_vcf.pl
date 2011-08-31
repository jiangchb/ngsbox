#!/usr/bin/perl

use strict;
use warnings;

my $vcf = shift or die "Please specify VCF file\n";

open VCF, $vcf or die "Cannot open input file\n";

# Loop through original reads and get prb entry
while( <VCF> ) {
	chomp;

	if($_ =~ /^#/) {
		print "$_\n";
	}
	else {
		my @a = split("\t", $_);
		my @b = split(";", $a[7]);

		my ($junk, $freq) = split("=", $b[2]);

		$a[7] =~ s/AC/DP/g;

		if($a[2] =~ /ins/) {
			$a[1]--;
		}

		if($freq > 0.8) {
			print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t" . "GT\t" . "1/1\n";
		}
		else {
			print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t" . "GT\t" . "0/1\n";
		}
	}
}

close VCF;

exit(0);
