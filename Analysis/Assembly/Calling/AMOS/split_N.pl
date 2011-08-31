#!/usr/bin/perl

use strict;
use warnings;

### User params
my $contig_fasta = shift;
open CONTIG, $contig_fasta or die "Cannot open $contig_fasta\n";

my $counter = 1;
while(<CONTIG>) {
	chomp($_);

	if (substr($_, 0, 1) eq ">") {
		my $id = substr($_, 1);
		my $seq = <CONTIG>;
		chomp($seq);

		if( ($seq =~ /NNN/) && ($id =~ /velvet/) ) {
			my @sub_seqs = split(/N+/, $seq);
			foreach my $sub_seq (@sub_seqs) {
				print ">$counter$id\n$sub_seq\n";
				$counter++;
			}
		}
		else {
			print ">$id\n$seq\n";
		}
	}
}
close CONTIG;
