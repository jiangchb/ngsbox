#!/usr/bin/perl
use strict;
use warnings;
use DBI;

my $usage = "\n$0 snp_calls ref_fasta\n\n";
my $snp   = shift or die $usage;
my $fasta = shift or die $usage;

my $refseq  = ();
my %refcall = ();


### Read fasta sequence
open FASTA, $fasta or die "Cannot open $fasta\n";
while(<FASTA>) {
	chomp;
	if($_ !~ />/) { 
		$refseq .= $_; 
	}
}

open SNP, $snp or die "Cannot open $snp\n";
my $target_start = 450;
my $target_length = 100;

while(<SNP>) {
        chomp;
	my @a = split(/\t/, $_);
	my $validation_start = $a[2] - 500;

	# Use original reference sequence for primer design
	print ">$a[0]-$a[1]-$a[2]-$a[3]($a[4]):TARGET=$target_start,$target_length\n";
	print substr($refseq, $validation_start, 800)  . "\n";
}
