#! /usr/bin/perl

use strict;

my $scaffile = shift or die "\n$0 scaffolds infolder chr\n\n";
my $infolder = shift or die "\n$0 scaffolds infolder chr\n\n";
my $chr = shift or die "\n$0 scaffolds infolder chr\n\n";

### Parse scaffolds
open SCAFF, $scaffile or die "\n$0 scaffolds infolder chr\n\n";
my %scaff_length = ();
while(<SCAFF>) {
	chomp;
	if( $_ =~ />/ ) {
		my $header = $_;
		my $seq = <SCAFF>;
		my ($scaff_id, $junk) = split(" | ", substr($header, 1));
		$scaff_length{$scaff_id} = length($seq);
		#print STDERR "$scaff_id\t" . $scaff_length{$scaff_id} . "\n";
	}

}

my @ali_scaff = glob("$infolder/alignments_Scaffold_*");

foreach my $infile ( @ali_scaff ) {
	
	open IN, $infile or die "Cannot open $infile\n\n";

	### Get scaffold name
	my $scaff_id = "NA";
	while(<IN>) {
		chomp;

		if($_ =~ /Alignments/) {
			my @e = split(" ", $_);
			$scaff_id = $e[5];
			last;
		}
	}

	### Find start
	while(<IN>) {
		chomp;
		
		if($_ =~ /BEGIN/) {
			<IN>;
			<IN>;
			my $refline = <IN>;
			my $scaffline = <IN>;
			my ($ref_pos, $refseq) = split(" ", $refline);
			my ($scaff_pos, $scaffseq) = split(" ", $scaffline);

			if($scaff_pos <= 75) {
				print "BEG\t$chr\t$ref_pos\t$scaff_pos\n";
			}

			last;
		}
	}

	### Find end
	my $ref_end    = -1;
	my $scaff_end = -1;
	my $junk = "";

	while(<IN>) {
		chomp;

		if($_ =~ /END/) {
			my @e = split(" - ", $_);
			my @r = split(" ", $e[1]);
			$ref_end = $r[0];
			($scaff_end, $junk) = split(" ", $e[2]);
		}
	}

	my $dist = $scaff_length{$scaff_id} - $scaff_end;

	if( ($scaff_id ne "NA")  && ($dist <= 75) ) {
		print "END\t$chr\t$ref_end\t$dist\n";
	}
}


