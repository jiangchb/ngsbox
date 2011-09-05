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
#  Module: Analysis::Eighties::Recomb_at_SV::calc_LD_diff.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $region_size = 1500;

my $max_del = 0;
my $min_called = 76;
my $min_allele_freq = 2;

my $database = "ath_eighty";

my $dbh;
&connect_to_db();

my $usage = "$0 PairedSiteFile\n";

my $file = shift or die $usage;

my $debug = 0;

########################################################################

open FILE, $file or die "Cannot open file: ".$file."\n";

while (<FILE>) {

	my @a = split " ";

	my $sv_chr = $a[0];
	my $sv_begin = $a[1];
	my $sv_end = $a[2];
	my $snp_chr = $a[3];
	my $snp_pos = $a[4];
	
	#print "##############################################################\n";
	#print $sv_chr, "\t", $sv_begin, "\t", $sv_end, "\t\t", $snp_chr, "\t", $snp_pos, "\t", $snp_pos+$region_size+1, "\n";

	# compare
	my $sv_avg_r2 = get_avg_r2_btw_regions($sv_chr, $sv_begin - $region_size - 1, $sv_begin - 1, $sv_chr, $sv_end + 1, $sv_end + $region_size + 1);
	my $snp_avg_r2 = get_avg_r2_btw_regions($snp_chr, $snp_pos - $region_size - 1, $snp_pos - 1, $snp_chr, $snp_pos + 1, $snp_pos + $region_size + 1);
	
	if ($sv_avg_r2 != -1 and $snp_avg_r2 != -1) {	
		print $sv_avg_r2, "\t", $snp_avg_r2, "\n";
	}

}

close FILE;

sub get_avg_r2_btw_regions {

	my ($chr1, $start1, $end1, $chr2, $start2, $end2) = @_; 

	my $region1_arr;
	my $region2_arr;

	$region1_arr = get_snps($chr1, $start1, $end1);
	$region2_arr = get_snps($chr2, $start2, $end2);
	
	
	#print "Num SNPs (left, right): ", @{$region1_arr}+0, "\t", @{$region2_arr}+0, "\n";

	my $r2_sum = 0;
	my $r2_count = 0;

	foreach my $snppos1 (@{$region1_arr}) {

		my $chr1 = shift @{$snppos1};
	        my $pos1 = shift @{$snppos1};
        	my $alleles1 = shift @{$snppos1};

		foreach my $snppos2 (@{$region2_arr}) {
			
			my $chr2 = shift @{$snppos2};
		        my $pos2 = shift @{$snppos2};
        		my $alleles2 = shift @{$snppos2};

			#print join("", @{$snppos1}), "<-----1\n";
			#print join("", @{$snppos2}), "<-----2\n";
			#print $alleles1, "\t", $alleles2, "\n";
			my $r2 = get_r2($alleles1, $alleles2, $snppos1, $snppos2);
			#print $r2, "\n";
			
			if ($r2 != -1) {
				$r2_sum += $r2;
				$r2_count++;
			}
		
			unshift @{$snppos2}, $alleles2;
			unshift @{$snppos2}, $pos2;
			unshift @{$snppos2}, $chr2;
		}
	}

	if ($r2_count != 0) {
		return $r2_sum/$r2_count;
	}
	else {
		return -1;
	}
	
}

sub get_r2 {
	my ($alleles1, $alleles2, $arr1, $arr2) = @_;

	if (length($alleles1) == 2 and length($alleles2) == 2) {
		my ($x11, $x12, $x21, $x22) = get_haplotype_freq($alleles1, $alleles2, $arr1, $arr2);

		my $D = $x11 - (($x11 + $x12) * ($x11 + $x21));

		#print substr($alleles1, 0, 1), substr($alleles1, 0, 1), ":", $x11, "\t", substr($alleles1, 0, 1), substr($alleles1, 1, 1), ":", $x12, "\t";
		#print substr($alleles1, 1, 1), substr($alleles2, 0, 1), ":", $x21, "\t", substr($alleles1, 1, 1), substr($alleles2, 1, 1), ":", $x22, "\n";

		if (($x11 + $x12) != 0 and ($x21 + $x22) != 0 and ($x11 + $x21) != 0 and ($x12 + $x22) != 0) {
			my $r = $D / sqrt(($x11 + $x12)*($x21 + $x22)*($x11 + $x21)*($x12 + $x22));
			return $r*$r;
		}
	}
	return -1;	

}

sub get_haplotype_freq {
	my ($alleles1, $alleles2, $arr1, $arr2) = @_;

	my $A1 = substr($alleles1, 0, 1);
	my $A2 = substr($alleles1, 1, 1);
	my $B1 = substr($alleles2, 0, 1);
	my $B2 = substr($alleles2, 1, 1);
		
	my $x11 = 0;
	my $x12 = 0;
	my $x21 = 0;
	my $x22 = 0;
	my $sum = 0;

	for (my $i = 0; $i < @{$arr1}; $i++) {
		if (${$arr1}[$i] eq $A1) {
			if (${$arr2}[$i] eq $B1) {
				$x11++;
				$sum++;
			}
			elsif (${$arr2}[$i] eq $B2) {
				$x12++;
				$sum++;
			}
		}
		elsif (${$arr1}[$i] eq $A2) {
			if (${$arr2}[$i] eq $B1) {
				$x21++;
				$sum++;
			}
			elsif (${$arr2}[$i] eq $B2) {
				$x22++;
				$sum++;
			}
		}
	}
	
	return ($x11/$sum, $x12/$sum, $x21/$sum, $x22/$sum);

}



sub get_snps {
	my ($chr, $begin, $end) = @_;	
	
	my $q = "
SELECT chromosome, position, ref, Agu_1_call as Agu_1, Bak_2_call as Bak_2, Bak_7_call as Bak_7, Cdm_0_call as Cdm_0, Del_10_call as Del_10, Dog_4_call as Dog_4, Don_0_call as Don_0, Ey15_2_call as Ey15_2, Fei_0_call as Fei_0, HKT2_4_call as HKT2_4, ICE102_call as ICE102, ICE104_call as ICE104, ICE106_call as ICE106, ICE107_call as ICE107, ICE111_call as ICE111, ICE112_call as ICE112, ICE119_call as ICE119, ICE120_call as ICE120, ICE127_call as ICE127, ICE130_call as ICE130, ICE134_call as ICE134, ICE138_call as ICE138, ICE150_call as ICE150, ICE152_call as ICE152, ICE153_call as ICE153, ICE163_call as ICE163, ICE169_call as ICE169, ICE173_call as ICE173, ICE181_call as ICE181, ICE1_call as ICE1, ICE212_call as ICE212, ICE213_call as ICE213, ICE216_call as ICE216, ICE21_call as ICE21, ICE226_call as ICE226, ICE228_call as ICE228, ICE29_call as ICE29, ICE33_call as ICE33, ICE36_call as ICE36, ICE49_call as ICE49, ICE50_call as ICE50, ICE60_call as ICE60, ICE61_call as ICE61, ICE63_call as ICE63, ICE70_call as ICE70, ICE71_call as ICE71, ICE72_call as ICE72, ICE73_call as ICE73, ICE75_call as ICE75, ICE79_call as ICE79, ICE7_call as ICE7, ICE91_call as ICE91, ICE92_call as ICE92, ICE93_call as ICE93, ICE97_call as ICE97, ICE98_call as ICE98, Istisu_1_call as Istisu_1, Kastel_1_call as Kastel_1, Koch_1_call as Koch_1, Lag2_2_call as Lag2_2, Leo_1_call as Leo_1, Lerik1_3_call as Lerik1_3, Memrut_1_call as Memrut_1, Mer_6_call as Mer_6, Nie1_2_call as Nie1_2, Ped_0_call as Ped_0, Pra_6_call as Pra_6, Qui_0_call as Qui_0, Rue3_1_call as Rue3_1, Sha_call as Sha, Star_8_call as Star_8, TueSB30_call as TueSB30, Tuescha9_call as Tuescha9, TueV12_call as TueV12, TueWa1_2_call as TueWa1_2, Vash_1_call as Vash_1, Vie_0_call as Vie_0, WalhaesB4_call as WalhaesB4, Xan_1_call as Xan_1, Yeg_1_call as Yeg_1 
FROM genome_matrix_calls_25_10_final_idx 
WHERE chromosome = $chr and position between $begin and $end 
ORDER BY chromosome, position";

        my $sth = $dbh->prepare($q) or die "cannot prepare\n";
        $sth->execute() or die "cannot execute\n";

	my @snps = ();

	while(my $ref = $sth->fetchrow_hashref()) {
		my @snp = ();

		my $count_a = 0;
		my $count_c = 0;
		my $count_g = 0;
		my $count_t = 0;
		my $count_del = 0;

		my $allele;

		$allele = $ref->{Agu_1};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Bak_2};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Bak_7};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Cdm_0};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Del_10};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Dog_4};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Don_0};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Ey15_2};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Fei_0};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{HKT2_4};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE102};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE104};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE106};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE107};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE111};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE112};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE119};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE120};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE127};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE130};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE134};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE138};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE150};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE152};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE153};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE163};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE169};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE173};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE181};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE1};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE212};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE213};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE216};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE21};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE226};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE228};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE29};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE33};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE36};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE49};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE50};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE60};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE61};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE63};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE70};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE71};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE72};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE73};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE75};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE79};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE7};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE91};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE92};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE93};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE97};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{ICE98};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Istisu_1};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Kastel_1};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Koch_1};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Lag2_2};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Leo_1};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Lerik1_3};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Memrut_1};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Mer_6};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Nie1_2};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Ped_0};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Pra_6};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Qui_0};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Rue3_1};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Sha};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Star_8};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{TueSB30};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Tuescha9};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{TueV12};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{TueWa1_2};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Vash_1};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Vie_0};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{WalhaesB4};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Xan_1};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";
		$allele = $ref->{Yeg_1};
		push @snp, $allele;
		$count_del++ if $allele eq "D" or $allele eq "-";
		$count_a++ if $allele eq "A";
		$count_c++ if $allele eq "C";
		$count_g++ if $allele eq "G";
		$count_t++ if $allele eq "T";

		if (	$count_del == $max_del and
			($count_a + $count_c + $count_g + $count_t) >= $min_called
			)
		{
			my $alleles = "";
			if ($count_a >= $min_allele_freq) { 
				$alleles .= "A";
			}
			if ($count_c >= $min_allele_freq) {
				$alleles .= "C";
			}
			if ($count_g >= $min_allele_freq) {
				$alleles .= "G";
			}
			if ($count_t >= $min_allele_freq) {
				$alleles .= "T";
			}
	
			if (length($alleles) >= 2) {	
				unshift @snp, $alleles;
		                unshift @snp, $ref->{position};
				unshift @snp, $ref->{chromosome};
				push @snps, \@snp;
			}
		}
	}


	$sth->finish();

	return \@snps;

}







sub connect_to_db {
        my $databaseName = $database;
        my $driver = "mysql";
        my $host = "orb.eb.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) || db_connect_error();
}


