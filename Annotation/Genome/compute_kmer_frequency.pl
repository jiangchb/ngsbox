#! /usr/bin/perl
use strict;
use warnings;

###### 
# NGSbox - bioinformatics analysis tools for next generation sequencing data
#
# Copyright 2007-2011 Stephan Ossowski, Korbinian Schneeberger
# 
# NGSbox is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or any later version.
#
# NGSbox is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# Please find the GNU General Public License at <http://www.gnu.org/licenses/>.
#
#  -------------------------------------------------------------------------
#
#  Module: Annotation::Genome::compute_kmer_frequency.pl
#  Purpose:
#  In:
#  Out:
#


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
