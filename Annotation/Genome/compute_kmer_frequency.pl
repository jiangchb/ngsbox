#!/usr/bin/perl
use strict;
use warnings;

my $usage = "\n\n$0 kmer_min_length kmer_max_length fasta_file\n\n";

my $kmer_min = shift or die $usage;
my $kmer_max = shift or die $usage;
my $file     = shift or die $usage;

open FILE, $file or die "Cannot open input file\n";;

### Count kmers
my %kmers = ();
my $glen = 0;

while( <FILE> ) {
	chomp;
	if($_ !~ />/) {
		my $seq = $_;
		for(my $i = 0; $i < length($seq) - $kmer_max; $i++) {
			for(my $kmer_len = $kmer_min; $kmer_len <= $kmer_max; $kmer_len++) {
				my $kmerseq = substr($seq, $i, $kmer_len);
				$kmers{$kmerseq}++;
			}
		}
		$glen += length($seq);
	}
}

close FILE; 


### Print kmer frequency
foreach my $kmerseq (sort keys %kmers) {
	if($kmerseq !~ /N/) {
		my $freq = $kmers{$kmerseq} / $glen;
		print $kmerseq ."\t". $kmers{$kmerseq} ."\t". $freq ."\n";
	}
}

exit(0);
