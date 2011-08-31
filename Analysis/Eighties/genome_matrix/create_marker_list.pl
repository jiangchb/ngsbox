#! /usr/bin/perl
use strict;
use warnings;

# Idea is to select position which are (almost) unique for the accession
# for which the quality_variant file is provided. These are good candidates
# for being markers for a mapping with an unsequenced accession.

my $usage = "\n$0 genome_matrix.frequency.txt quality_variant.txt quality_threshold\n\n";
my $freq = shift or die $usage;
my $variants = shift or die $usage;
my $qthres = shift or die $usage;

open VAR, $variants or die $usage;

my %MARKER = ();
my $sample = "";
while (my $line = <VAR>) {
	my @a = split " ", $line;
	$sample = $a[0];
	if ($a[4] ne "-" and $a[5] >= $qthres and $a[8] == 1) {
		my $nuc = 0;
		if ($a[4] eq "A") {
			$nuc = 0;
		}
		elsif ($a[4] eq "C") {
			$nuc = 1;
                }
		elsif ($a[4] eq "G") {
			$nuc = 2;
                }
		elsif ($a[4] eq "T") {
			$nuc = 3
                }
		else {
			print STDERR $line;
			exit(1);
		}
		$MARKER{$a[1]."#".$a[2]} = $nuc;
	}
}

close VAR;

open FREQ, $freq or die $usage;
my $c = 0;
while (my $line = <FREQ>) {
	my @a = split " ", $line;
	
	print STDERR $a[0], "\t", $a[1], "\n" if ($a[1]%100000) == 0;
	
	if (defined($MARKER{$a[0]."#".$a[1]})) {
		
		my $count_n = $a[13];
		if ($count_n <= 60) {  # minimum 60 genomes called
			if ($a[8+$MARKER{$a[0]."#".$a[1]}] != 0) { # accession allele was seen in the 80 genomes
				delete $MARKER{$a[0]."#".$a[1]};
			}
			else {
				my $count_alleles = 0;
				$count_alleles++ if $a[8] > 0;
				$count_alleles++ if $a[9] > 0;
				$count_alleles++ if $a[10] > 0;
				$count_alleles++ if $a[11] > 0;

				if ($count_alleles > 1) {
					delete $MARKER{$a[0]."#".$a[1]};
				}
			}
		}
	}
}

open OUT, "> putative_markers.$sample.$qthres.tmp";
foreach my $key (sort keys %MARKER) {
	my @a = split "#", $key;
	print OUT $sample, "\t", $a[0], "\t", $a[1],"\n";
}
close OUT;

system("sort -n -k2 -k3 putative_markers.$sample.$qthres.tmp > putative_markers.$sample.$qthres.txt");
system("rm putative_markers.$sample.$qthres.tmp");




