#!/usr/bin/perl 

####################################################################
# Overlap Genotyping SNPs with NGS SNPs
# Author: Stephan Ossowski
####################################################################

use strict;
use warnings;

my $usage = "\n\n$0 enrichment_file SNP_VCF\n\n";
my $enriched = shift or die $usage;
my $snp_vcf  = shift or die $usage;


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



# Read VCF and return SNP from enriched regions
open VCF, $snp_vcf or die "Cannot open $snp_vcf file\n";
while( <VCF> ) {

	chomp;

	if($_ =~ /#/) {
		print "$_\n";
	}
	else {
		my @a = split("\t", $_);

		if( ! exists $enr{$a[0]}{$a[1]} ) {
			$a[6] .= ",NOTENRICHED";
		}

		print join("\t",@a);
		print "\n";
		#for (my $i = 0; $i < $#a; $i++) {
		#	print "$a[$i]\t";
		#}
	}
}

close VCF;
exit(0);
