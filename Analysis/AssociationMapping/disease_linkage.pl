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
#  Module: Analysis::AssociationMapping::disease_linkage.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "\n$0 win_size win_step vcf_files<comma separated\n\n";

my $win_size  = shift or die $usage;
my $win_step  = shift or die $usage;
my $vcf_files = shift or die $usage;

my @vcf_file = split(",", $vcf_files);

my %snps = ();
my %snp_summary = ();
my $sample_counter = 0;

foreach my $vcf ( @vcf_file ) {

	$sample_counter++;

	open VCF, $vcf or die "Cannot open input file\n";

	while( my $line = <VCF> ) {
		chomp($line);

		if($line =~ /^#/) {
			# ignore for now
		}
		else {
			my @a        = split("\t", $line);
			my @info     = split(";", $a[7]);	# Generic info about call
			my @gt_order = split(":", $a[8]);	# Genotypes order
			my @gt_value = split(":", $a[9]);	# Genotypes values

			# Get infotation
			my %info = ();
			foreach my $info_string (@info) {
				my($type, $value) = split("=", $info_string);
				$info{$type} = $value;
			}


			### Get genotype
			my %geno = ();
			for(my $i = 0; $i <= $#gt_order; $i++) {
				$geno{$gt_order[$i]} = $gt_value[$i];
			}


			### Allele support
			my $quality_cov = "\\N";
			my $ref_support = "\\N";
			my $snp_support = "\\N";
			if(exists $geno{DP}) {
				$quality_cov = $geno{DP}
			}
#			if(exists $geno{AD}) {
#				($ref_support, $snp_support) = split(",", $geno{AD});
#			}


			### Enriched regions
			my $enriched = 1;
			if($a[6] =~ /NOTENRICHED/) {
				$enriched = 0;
			}


			### Filter
			my $passed_filter = 0;
			if($a[6] =~ /PASS/) {
				$passed_filter = 1;
			}


			### Hom/Het
			my $homhet = "het";
			if($geno{GT} eq "1/1") {
				$homhet = "hom";
			}


			### Store SNP
			$snps{$sample_counter}{$a[0]}{$a[1]} = $homhet;

			if( ($info{set} eq "Intersection") && ($homhet eq "hom") && ($quality_cov >= 15) ) {
				$snp_summary{$a[0]}{$a[1]}++;
			}
		}
	}

	close VCF;
}
print STDERR "\n\nFinished reading VCF files\n\n";



### Compute shared SNPs per region
foreach my $chr ( sort {$a<=>$b} keys %snp_summary ) {

	### Compute max SNP
	my $max_size = 0;
	foreach my $pos ( keys %{$snp_summary{$chr}} ) {
		if($pos > $max_size) {
			$max_size = $pos;
		}
	}

	### Move sliding window along chromosome
	for(my $win_start = 1; $win_start < $max_size; $win_start += $win_step) {

		my $win_end = $win_start + $win_size - 1;
		my $total_shared = 0;
		my $total_tests = 0;
		
		foreach my $pos ( keys %{$snp_summary{$chr}} ) {

			if( ($pos >= $win_start) && ($pos <= $win_end) && ($snp_summary{$chr}{$pos} < $sample_counter) ) {
				my $occurrence = 0;
				my $homness = 0;
				
				for( my $sample = 1; $sample <= $sample_counter; $sample++ ) {
					
					if( exists $snps{$sample}{$chr}{$pos} ) {
						$occurrence++;

						if($snps{$sample}{$chr}{$pos} eq "hom") {
							$homness = 1;
						}
					}
				}

				if($homness == 1) {
					$total_tests += $sample_counter;
					$total_shared += $occurrence;
				}
			}
		}

		if( ($total_tests / $sample_counter) >= ($win_size / 1000000 * 5) ) {
			my $freq = $total_shared / $total_tests;

			### Smooth frequency
			$freq += 0.01;
			if($freq > 1.0) { $freq = 1.0; }
			my $score = $freq ** 4;

			### Print score (smoothed frequency)
			if( $score > 0.95 ) {
				print "$chr\t$win_start\t$win_end\t$total_shared\t$total_tests\t$score\n";
			}

			### For R plotting
			# my $tmp_pos = 1000000000 + $win_start;
			# if($chr eq "X") { $chr = 23; }
			# elsif($chr eq "Y") { $chr = 24; }
			# elsif($chr eq "MT") { $chr = 25; }
			# my $win_start_plot = $chr . $tmp_pos;
			# print "$win_start_plot\t$score\n";
		}
	}
}

exit(0);
