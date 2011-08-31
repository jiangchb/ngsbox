#!/usr/bin/perl
use strict;
use warnings;
use DBI;

my $flat  = shift;
my $fasta = shift;

my $id = "";
my %contigs= ();

open FASTA, $fasta or die "Cannot open $fasta\n";
while(<FASTA>) {
	chomp;
	if($_ =~ />/) { $id = $_; }
	else { $contigs{$id} = $_; }
}

open FLAT, $flat or die "Cannot open $flat\n";
while(<FLAT>) {
        chomp;
	my @elem = split(/\t/, $_);
	my @snps = split(/;/, $elem[4]);
	my $target_start;
	my $target_length;

	if(scalar(@snps) == 1) {
		$target_start = $snps[0] - $elem[2] - 50;
		$target_length = 100;
	}
	elsif(scalar(@snps) > 1) {
		$target_start = $snps[0] - $elem[2] - 50;
		$target_length = $snps[-1] - $snps[0] + 100;
		if( 	($target_start < 30) ||
			($target_length > 400) || 
			($target_length < 0) 
		) { print STDERR "ERROR:$elem[0]-$elem[2]-$elem[3]\n"; }
	}
	else { print STDERR "ERROR\n"; }

	print ">$elem[0]-$elem[1]-$elem[2]-$elem[3]($elem[4]):TARGET=$target_start,$target_length\n";
	print $contigs{">$elem[1]-$elem[2]-$elem[3]"} . "\n";
}
