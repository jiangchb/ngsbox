#!/usr/bin/perl
use strict;
use warnings;
use DBI;

my $usage = "\n$0 regions_file ref_fasta\n\n";
my $regions   = shift or die $usage;
my $fasta = shift or die $usage;

my %refseq  = ();
my $chr = "";

### Read fasta sequence
open FASTA, $fasta or die "Cannot open $fasta\n";
while(<FASTA>) {
	chomp;

	if($_ =~ />/) {
		$chr = substr($_, 1);
	}
	else { 
		$refseq{$chr} .= $_; 
	}
}


open REGION, $regions or die "Cannot open $regions\n";

while(<REGION>) {
        chomp;
	my @a = split(/\t/, $_);
	my $chr = $a[1];
	my $length = $a[3] - $a[2] + 1;
	my $middle = $a[2] + int($length / 2);
	my $validation_start = $middle - 75;
	my $validation_end = $validation_start + 150;
	

	print ">$chr-$validation_start-$validation_end($middle):TARGET=75,1\n";
	print substr($refseq{$chr}, $validation_start, 150)  . "\n";
}
