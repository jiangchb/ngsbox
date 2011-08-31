#! /usr/bin/perl

use lib "$ENV{PGSP}/Prediction/Validation/";
use Genome;
use strict;

my $usage = "$0 chrfasta snpfile delfile insfile chr begin end\n";

my $reffile = shift or die $usage;
my $snpfile = shift or die $usage;
my $deletionfile = shift or die $usage;
my $insertionfile = shift or die $usage;
my $chr = shift or die $usage;
my $begin = shift or die $usage;
my $end = shift or die $usage;

my $genome = new Genome();
$genome->calc_chr_seq($reffile, $snpfile, $deletionfile, $insertionfile);

my $seq = $genome->get_sequence($chr, $begin, $end);

print ">", $chr, " ", $begin, " ", $end, "\n";
for (my $i = 1; $i<=length($seq); $i++) {
	print substr($seq, $i-1, 1);
	print "\n" if $i%60 == 0;
}




