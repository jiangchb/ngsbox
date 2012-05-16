#! /usr/bin/perl

use strict;

my $usage = "$0  variation-file  reference-file  mutation/marker-file\n";
my $varfile = shift or die $usage;
my $reffile = shift or die $usage;
my $mutfile = shift or die $usage;

my $MIN_QUAL_CALLED = 40;
my $MAX_REPETITIVITY_CALLED = 1.0;
my $REGION_SIZE = 30;

my $NUM_NOT_CALLED = 1;

my %PASSED_FILTER = ();

parse($varfile);
parse($reffile);


open FILE, $mutfile or die "cannot open file $mutfile\n";
my $outfile = $mutfile;
$outfile =~ s/\.txt//g;
open OUT, ">$outfile.VICINITY.q$MIN_QUAL_CALLED.rep$MAX_REPETITIVITY_CALLED.size$REGION_SIZE.match$NUM_NOT_CALLED.txt" or die "cannot open output file\n";

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[5] >= $MIN_QUAL_CALLED and $a[8] <= $MAX_REPETITIVITY_CALLED) {
		my $chr = $a[1];
		my $pos = $a[2];
		my $count = 0;
		for (my $i = $pos - $REGION_SIZE; $i <= $pos + $REGION_SIZE; $i++) {
			$count++ if not defined($PASSED_FILTER{$chr}{$i});
		}
		if ($count <= $NUM_NOT_CALLED) {
			print OUT $line;
		}
	}
}




sub parse {
	my ($file) = @_;

	my $c = 0;
	open FILE, $file or die "cannot open file $file\n";
	while (<FILE>) {
		$c++;
		print STDERR $_ if $c%400000==0;
		my @a = split " ";
		if ($a[5] >= $MIN_QUAL_CALLED and $a[8] <= $MAX_REPETITIVITY_CALLED) {
			$PASSED_FILTER{$a[1]}{$a[2]} = 1;
		}
	}
	close FILE;
}


