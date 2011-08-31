#! /usr/bin/perl
use strict;
use warnings;

my $usage = "$0 <fastafile> <outfile prefix>\n";

my $in  = shift or die $usage;
my $out = shift or die $usage;

open IN, $in or die $usage;
open OUT1, ">$out.1" or die $usage;
open OUT2, ">$out.2" or die $usage;


while (my $header1 = <IN>) {

	my $seq1    = <IN>;
	my $header2 = <IN>;
	my $seq2    = <IN>;

	chomp($header1);
	chomp($header2);
	chomp($seq1);
	chomp($seq2);


	# Quality string
	my $q2 = "";

	for (my $i = 0; $i < length($seq1); $i++) {
		$q2 .= "I";
	}

	print OUT1 substr($header1, 1) . "\t$seq1\t1\t$q2\n";
	print OUT2 substr($header1, 1) . "\t$seq2\t2\t$q2\n";
}
