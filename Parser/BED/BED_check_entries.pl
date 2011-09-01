#!/usr/bin/perl
use strict;
use warnings;

my $file = shift;

my $last_chr = "";
my %cds = ();

open IN, $file or die "Cannot open input file\n";


while( <IN> ) {
	chomp;
	my ($chr, $beg, $end, $strand) = split(/\t/, $_);

	# New chromosome: reset CDS hash
	if($chr ne $last_chr) {
		$last_chr = $chr;
		%cds = ();
	}

	if($beg <= $end) {
		if( ! exists $cds{$beg} ) {
			print "$_\n";
			$cds{$beg} = $end;
		}
		elsif( $cds{$beg} != $end ) {
			print "$_\n";
			$cds{$beg} = $end;
		}
	}
	else {
		print STDERR "$_\n";
	}
}

close IN;

exit(0);
