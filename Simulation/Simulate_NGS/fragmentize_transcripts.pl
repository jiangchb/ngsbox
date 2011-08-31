#!/usr/bin/perl 
##############################################
# Simulate Illumina GA2/HiSeq RNA-seq reads
# DEPRECATED! Use METASIM instead.
##############################################

use strict;
use warnings;

my $cDNA_file       = shift;
my $expression_file = shift;
my $frag_length     = shift;
my $num_reads       = shift;

### Read cDNA file in fasta format
open(TRANS, $cDNA_file ) or die "Cannot open input file $cDNA_file\n";
my $seq;
my %RNA = ();
while( <TRANS> ) { 
	chomp; 
	if(substr($_, 0, 1) ne ">") { 
		$seq .= $_; 
	} else {
		my @a = split " ";
		$RNA{substr($a[0], 1, length($a[0])-1)} = $seq;
		$seq = ""; 
	}
}
close TRANS;


### Read expression profile
open(EXPRESS, $expression_file ) or die "Cannot open input file $expression_file\n";
my %EXPR = ();
my $complete_expr = 0;
while(<EXPRESS>) {
	my @a = split " ";
	my $id = $a[0].".".$a[1];
	$EXPR{$id} = $a[3];
	$complete_expr += $a[3];
}
close EXPRESS;


### Join cDNA and expression profiles
foreach my $transcript (keys %RNA) {
	if (defined($EXPR{$transcript})) {
		$EXPR{$transcript} = $EXPR{$transcript} / $complete_expr;
	} else {
		$EXPR{$transcript} = 0;
	}
}


### Create random set of Illumina reads from transcripts
my $count = 100000000;
foreach my $transcript (keys %RNA) {

	my $expr = $EXPR{$transcript};
	my $seq = $RNA{$transcript};

	# Randomize fragment
	if (defined($seq) and length($seq) >= $frag_length) {	
		my $num_read_frag = int($num_reads * $expr);
		
		for (my $i = 0; $i<$num_read_frag; $i++) {
			my $start = int(rand(length($seq) - $frag_length));
			my $read = substr($seq, $start, $frag_length);
	
			for(my $i = 0; $i < $frag_length; $i++) {
				my $base = substr($read, $i, 1);
				$base = &rand_base($base);
				substr($read, $i, 1, $base);
			}

			if(int(rand(2)) == 1) {
				my $rev_seq = reverse($read);
				$rev_seq =~ tr/acgtACGT/TGCATGCA/;
				$read = $rev_seq;
			}
		
			if($read =~ /^[ATCG]+$/) {
				print "f$count\t$read\n";
			}
			$count++;
		}
	}
}

exit(0);


### Random sequencing error
sub rand_base
{
	my $base = shift;
	my @nucleotides = ('A', 'T', 'C', 'G');
	my $error_base = $base;
	my $error_rate = 100;
	my $error_event = int(rand($error_rate));
	if($error_event < 2) {
		while($error_base eq $base) {
			my $base_event = int(rand(4));
			$error_base = $nucleotides[$base_event];
		}
	}
	return($error_base);
}
