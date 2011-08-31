#! /usr/bin/perl
use strict;
my $usage = "$0 gff [gene|transposable_element_gene]\n";

my $file = shift or die $usage;
my $feature = shift or die $usage;
open FILE, $file or die "Cannot open file\n";

print "#?chr\tpos\tsize\tstrand\tid\n";

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[2] eq $feature) {
		my @b = split "=", $a[8];
		my $id = substr($b[1], 0, 9);
		my $chr = substr($a[0], 3, 1);
		if ($chr eq "1" or $chr eq "2" or $chr eq "3" or $chr eq "4" or $chr eq "5") {
			print $chr, "\t", $a[3], "\t", $a[4]-$a[3]+1, "\t.\t", $id, "\n";
		}
	}
}


