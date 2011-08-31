#!/usr/bin/perl 

#############################################
# Simulate Illumina GA2/HiSeq DNA-seq reads
# DEPRECATED! Use METASIM instead.
#############################################

use strict;
use warnings;

my $usage = "perl fragmentize.pl <error flag 1|0> <multiple fasta file> <read length> <read number>";
my $with_errors_flag = shift or die $usage;
my $file = shift or die $usage;
my $frag_length = shift or die $usage;
my $frag_number = shift or die $usage;

my @seqs = ();
my $count = 100000000;

open FILE, $file or die $usage;

### Read sequences
my $id = "";
my $seq = "";
while(<FILE>) {
	chomp();
	if (substr($_, 0, 1) eq ">") {
		if ($id ne "") {
			push @seqs, $seq;
		}
		my @a = split " ";
		$id = substr($a[0], 1, length($a[0])-1);
		$seq = "";
	} else {
		if ($id ne "") {
			$seq .= $_;
		}
	}	
}
if ($id ne "") {
	push @seqs, $seq;
}

### Create random set of sequenced fragments from pseudo-genome
for( my $i = 0; $i < $frag_number; $i++ ) {
	
	# randomize contig
	my $contig = int(rand($#seqs+1));

	# cut out fragment
	my $start = int(rand(length($seqs[$contig]) - $frag_length - 2));
	my $seq = substr($seqs[$contig], $start, $frag_length);

	# introduce errors
	my $frag_seq = $seq;
	if ($with_errors_flag == 1) {
		$frag_seq = '';	
		for(my $i = 0; $i < $frag_length; $i++) {
			my $base = substr($seq, $i, 1);
			$base = &rand_base($base, $i);
			$frag_seq .= $base;
		}
	}

	# reverse every second fragment
	if( int(rand(2)) == 1) {
		my $rev_seq = reverse($frag_seq);
		$rev_seq =~ tr/acgtACGT/TGCATGCA/;
		$frag_seq = $rev_seq;
	}
	
	if($frag_seq =~ /^[ATCG]+$/) {
		print "f$count\t$frag_seq\n";
	}
	$count++;
}

exit(0);

### Random sequencing error
sub rand_base
{
	my $base = shift;
	my $position = shift;
	my @nucleotides = ('A', 'T', 'C', 'G');
	my $error_base = $base;
	my $error_rate = 100;
	if($position > 28) { $error_rate = 50; }
	if($position > 34) { $error_rate = 30; }
	my $error_event = int(rand($error_rate));
	if($error_event == 25) {
		while($error_base eq $base) {
			my $base_event = int(rand(4));
			$error_base = $nucleotides[$base_event];
		}
	}
	return($error_base);
}
