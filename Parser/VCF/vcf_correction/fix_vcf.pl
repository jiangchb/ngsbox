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
		my @b = split(":", $a[8]);
		my @c = split(":", $a[9]);

		print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t" . "$b[1]:$b[2]:$b[0]\t" . "$c[1]:$c[2]:$c[0]\n";
	}
}

close VCF;

exit(0);
