#! /usr/bin/perl
use strict;

my $usage ="$0 num_reads_file\n";
my $file = shift or die $usage;

my $MIN_NUM_READS = 200000;


open FILE, $file or die "cannot open file $file\n";

my %SIZE = ();

while (<FILE>) {
	my @a = split " ";
	$SIZE{$a[0]} = $a[1];
}

close FILE;



for (my $i = 1; $i <= 96; $i++) {

	if ($SIZE{$i} >= $MIN_NUM_READS) {
		chdir("sample_$i/map.list.corrected.sorted.Analysis/ConsensusAnalysis/supplementary_data");

		print("$i\n"); 
		system("perl /home/schneeberger/ngsbox/Analysis/GBS/gbs.sliding_window.pl consensus_summary.marker.txt /projects/dep_coupland/grp_nordstrom/sequencing/Athal/genome/strains/Ws-2-Tuebingen/AnalysisFolder_q01/ConsensusAnalysis/quality_variant.woDel.q32.cov6-39.txt /projects/dep_coupland/grp_nordstrom/data/Athal/TAIR10/chrsizes.chr1-5.txt");

		chdir("../../../../");
	}
}






