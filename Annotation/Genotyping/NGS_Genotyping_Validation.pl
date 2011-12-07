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
#  Module: Annotation::Genotyping::NGS_Genotyping_Validation.pl
#  Purpose:
#  In:
#  Out:
#


####################################################################
# Overlap Genotyping SNPs with NGS SNPs
# Author: Stephan Ossowski
####################################################################


my $usage = "\n\n$0 genotype_file enrichment_file, SNP_VCF\n\n";
my $genotype = shift or die $usage;
my $enriched = shift or die $usage;
my $snp_vcf  = shift or die $usage;

my %genotypes = ();


# Get enriched regions
my %enr = ();
open ENRICHED, $enriched or die "Cannot open $enriched file\n";
while( <ENRICHED> ) {
	chomp;
	my @a = split("\t", $_);
	if($a[3] eq "enriched") {
		for( my $i = $a[1]; $i <= $a[2]; $i++) {
			$enr{$a[0]}{$i} = 1;
		}
	}
}
close ENRICHED;


# Read NGS SNP file
my %snps   = ();
open SNPVCF, $snp_vcf or die "Cannot open $snp_vcf file\n";
while( <SNPVCF> ) {
	chomp;
	my @a = split("\t", $_);
	if(exists $enr{$a[0]}{$a[1]}) {
		$snps{$a[0]}{$a[1]} = $_;
	}
}
close SNPVCF;


# Get all genotyped positions and compare to NGS SNPs
my %stats = ();
open GENO, $genotype or die "Cannot open $genotype file\n";
while( <GENO> ) {
	chomp;
	my ($sample, $id, $chr, $pos, $allele_1, $allele_2, $refbase) = split("\t", $_);

	$allele_2 = substr($allele_2, 0, 1);

	### Filter genotypes (enriched region + clear call
	if( ($allele_1 =~ /[ACGT]/) && ($allele_2 =~ /[ACGT]/) && ($refbase =~ /[ACGT]/) && (exists $enr{$chr}{$pos}) ) {
		
		$stats{total_genotypes}++;

		### Case: NGS SNP called
		if( exists $snps{$chr}{$pos} ) {

			my @a   = split("\t", $snps{$chr}{$pos});
			my @AN  = split(";", $a[7]);      # Annotations
			my @GTo = split(":", $a[8]);      # Genotypes order
			my @GTv = split(":", $a[9]);      # Genotypes values

			$stats{total_ngsSNP}++;
			my $ngs_call = substr($a[4], 0, 1);


			# Case: Genotyping array calls SNP (True Positive)
			if( ($allele_1 ne $refbase) || ($allele_2 ne $refbase) ) {
				$stats{total_gSNP}++;
				$stats{TP}++;
			
				### Genotype array call is Het
				if($allele_1 ne $allele_2) {
					$stats{gSNP_het}++;
					$stats{TP_het}++;

					if( (($ngs_call eq $allele_1) || ($ngs_call eq $allele_2)) && ($GTv[0] eq "0/1") ) {
						$stats{good_het}++;
					}
					else {
						$stats{bad_het}++;
					}
				}
				### Genotype array call is Hom
				else {
					$stats{gSNP_hom}++;
					$stats{TP_hom}++;

					if( ($ngs_call eq $allele_1) && ($ngs_call eq $allele_2) && ($GTv[0] eq "1/1") ) {
						$stats{good_hom}++;
					}
					else {
						$stats{bad_hom}++;
					}
				}
			}
			# Case: Genotyping array calls Ref (False Positive)
			else {
				$stats{FP}++;

				if($GTv[0] eq "0/1") {
					$stats{FP_het}++;
				}
				elsif($GTv[0] eq "1/1") {
					$stats{FP_hom}++;
				}
			}
		}
		### Case: No NGS SNP called
		else {
			# Case: Genotyping array calls SNP (False negative)
			if( ($allele_1 ne $refbase) || ($allele_2 ne $refbase) ) {
				$stats{total_gSNP}++;
				$stats{FN}++;

				# Missed Het
				if($allele_1 ne $allele_2) {
					$stats{gSNP_het}++;
					$stats{FN_het}++;
				}
				# Missed Hom
				else {
					$stats{gSNP_hom}++;
					$stats{FN_hom}++;
				}
			}
			# Case: Genotyping array calls Ref (True Negative)
			else {
				$stats{TN}++;
			}
		}
	}
}

my $sensitivity = $stats{TP} / ($stats{TP} + $stats{FN});
my $sensitivity_het = $stats{TP_het} / ($stats{TP_het} + $stats{FN_het});
my $sensitivity_hom = $stats{TP_hom} / ($stats{TP_hom} + $stats{FN_hom});

my $PPV = $stats{TP} / ($stats{TP} + $stats{FP});
my $PPV_het = $stats{TP_het} / ($stats{TP_het} + $stats{FP_het});
my $PPV_hom = $stats{TP_hom} / ($stats{TP_hom} + $stats{FP_hom});


print "Sample: $genotype\n" .
	"\nTotal genotyped positions: " . $stats{total_genotypes} .
	"\nTotal genotyped SNPs: " . $stats{total_gSNP}++ .
	"\nTotal NGS SNPs: " . $stats{total_ngsSNP} .

	"\nNGS TP: " . $stats{TP} .
	"\nNGS FP: " . $stats{FP} .
	"\nNGS TN: " . $stats{TN} .
	"\nNGS FN: " . $stats{FN} .
	"\nNGS PPV: " . $PPV .
	"\nNGS Sensitivity: " . $sensitivity .
	"\n" .

	"\nGenotyping Het: " . $stats{gSNP_het} .
	"\nHet TP : " . $stats{TP_het} .
	"\nHet TP (correct allele): " . $stats{good_het} .
	"\nHet TP (wrong allele): " . $stats{bad_het} .
	"\nHet FP: " . $stats{FP_het} .
	"\nHet FN: " . $stats{FN_het} .
	"\nHet PPV: " . $PPV_het .
	"\nHet Sensitivity: " . $sensitivity_het .
	"\n".

	"\nGenotyping Hom: " . $stats{gSNP_hom} .
	"\nHom TP : " . $stats{TP_hom} .
	"\nHom TP (correct allele): " . $stats{good_hom} .
	"\nHom TP (wrong allele): " . $stats{bad_hom} .
	"\nHom FP: " . $stats{FP_hom} .
	"\nHom FN: " . $stats{FN_hom} .
	"\nHom PPV: " . $PPV_hom .
	"\nHom Sensitivity: " . $sensitivity_hom .
	"\n\n";

close GENO;
exit(0);
