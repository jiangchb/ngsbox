#!/usr/bin/perl 

####################################################################
# Overlap Genotyping SNPs with NGS SNPs
# Author: Stephan Ossowski
####################################################################

use strict;
use warnings;

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
				#print "$sample, $id, $chr, $pos, $allele_1, $allele_2\n" . $snps{$chr}{$pos} . "\n";
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
					$stats{missed_het}++;
				}
				# Missed Hom
				else {
					$stats{gSNP_hom}++;
					$stats{missed_hom}++;
				}
			}
			# Case: Genotyping array calls Ref (True Negative)
			else {
				$stats{TN}++;
			}
		}
	}
}
	
print "Sample: $genotype" .
	"\nTotal genotyped positions: " . $stats{total_genotypes} .
	"\nTotal genotyped SNPs: " . $stats{total_gSNP}++ .
	"\nTotal NGS SNPs: " . $stats{total_ngsSNP} .

	"\nNGS TP: " . $stats{TP} .
	"\nNGS FP: " . $stats{FP} .
	"\nNGS TN: " . $stats{TN} .
	"\nNGS FN: " . $stats{FN} .
	"\n" .
	"\nGenotyping Het: " . $stats{gSNP_het} .
	"\nMissed Het: " . $stats{missed_het} .
	"\nGood Het: " . $stats{good_het} .
	"\nBad Het: " . $stats{bad_het} .
	"\n".
	"\nGenotyping Hom: " . $stats{gSNP_hom} .
	"\nMissed Hom: " . $stats{missed_hom} .
	"\nGood Hom: " . $stats{good_hom} .
	"\nBad Hom: " . $stats{bad_hom} .
	"\n\n";

close GENO;
exit(0);
