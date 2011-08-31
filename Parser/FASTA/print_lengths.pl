#! /usr/bin/perl
use strict;
my $usage="$0 fasta\n";
open FILE, shift or die $usage;

my $seq = "";

while (<FILE>) {
	if (substr($_, 0, 1) eq ">") {
		print length($seq), "\n" if ($seq ne "");
		$seq = "";
	}
	else {
		chomp();
		$seq.=$_;
	}
}
print length($seq), "\n" if ($seq ne "");



