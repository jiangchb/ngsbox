#! /usr/bin/perl

use strict;

my $usage = "$0 variation-file reference-file mutation-file\n";
my $varfile = shift or die $usage;
my $reffile = shift or die $usage;
my $mutfile = shift or die $usage;

my $MIN_QUAL_CALLED = 25;
my $REGION_SIZE = 10;
my $NUM_NOT_CALLED = 0;

my %CALLED = ();

parse($varfile);
parse($reffile);


open FILE, $mutfile or die "cannot open file $mutfile\n";

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[5] >= $MIN_QUAL_CALLED) {
		my $chr = $a[1];
		my $pos = $a[2];
		my $count = 0;
		for (my $i = $pos - $REGION_SIZE; $i <= $pos + $REGION_SIZE; $i++) {
			$count++ if not defined($CALLED{$chr}{$pos});
		}
		if ($count <= $NUM_NOT_CALLED) {
			print $line;
		}
	}
}




sub parse {
	my ($file) = @_;

	open FILE, $file or die "cannot open file $file\n";
	while (<FILE>) {
		my @a = split " ";
		if ($a[5] >= $MIN_QUAL_CALLED) {
			$CALLED{$a[1]}{$a[0]} = 1;
		}
	}
	close FILE;
}


