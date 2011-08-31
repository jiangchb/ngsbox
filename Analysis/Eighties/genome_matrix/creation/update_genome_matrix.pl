#! /usr/bin/perl

my $usage = "\n$0 minQualityValue\n";

my $min_qual_value = shift or die $usage;

use strict;

my %DELETIONS = ();
my %UNSEQUENCED = ();


#chdir("/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2");
my $base = "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2";

my @ECOTYPES = ("Agu-1","Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169","ICE173","ICE181","ICE1","ICE212","ICE213","ICE216","ICE21","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE7","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2","Ped-0","Pra-6","Qui-0","Rue3-1","Sha","Star-8","TueSB30","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1");


for (my $i = 0; $i <@ECOTYPES; $i++) {
	my $ecotype = $ECOTYPES[$i];

	my $ecotype_folder = "/".$base."/".$ecotype."/";
	my $alignment_folder = $ecotype_folder."/AlignmentFolder/SV_lib1_version4_final/";
	my $sv_file = $alignment_folder."/SV_deletion_high_quality.lib1.PE.txt";

	# Read in SVdeletions
	print $sv_file, "\n";

	open SV, $sv_file or die "Cannot find file1: ",$sv_file, "\n";;
	my %DEL = ();
	while (my $line = <SV>) {
		my @a = split " ", $line;
		for (my $i = $a[4]; $i <= $a[5]; $i++) {
			$DEL{$a[3]."#".$i} = 1;
		}
	}
	close SV;

	$DELETIONS{$ecotype} = \%DEL;

	# Read in unseq
	my $folder = $ecotype_folder."/AlignmentFolder/Analysis_q01/ConsensusAnalysis";
	my $unseq_file = $folder."/unsequenced.txt";
	print  $unseq_file, "\n";

	open UNSEQ, $unseq_file or die "Cannot find file2: ", $unseq_file,"\n";
	my %ZERO = ();
        while (my $line = <UNSEQ>) {
                my @a = split " ", $line;
                for (my $i = $a[2]; $i <= $a[3]; $i++) {
                        $ZERO{$a[1]."#".$i} = 1;
                }
        }
        close UNSEQ;
	
	$UNSEQUENCED{$ecotype} = \%ZERO;

}


# For stats:
my %NUM_CALLS = ();
my %NUM_CALLS_WITH_SVDEL = ();

my $NUM_SNPS = 0;
my $NUM_SNPS_NONUNIQ = 0;

open GM, "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2/genome_matrix_new/genome_matrix_max_qual.25.$min_qual_value.txt"; 
open OUT, ">/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2/genome_matrix_new/genome_matrix_max_qual.25.$min_qual_value.wSVdel.wUnseq.txt";

while (my $line = <GM>) {
	my @a = split " ", $line;

	my $chr = $a[0];
	my $pos = $a[1];
	my $ref = $a[2];
	my $new_line = $chr."\t".$pos."\t".$ref;

	print STDERR $chr, "\t", $pos, "\n" if $pos%100000 == 0;
	
	my $num_calls = 0;
	my $num_calls_with_svdel = 0;
	my $num_snps = 0;

	for (my $e = 0; $e < @ECOTYPES; $e++) {
		if ($a[3+$e] eq "A" or $a[3+$e] eq "T" or $a[3+$e] eq "G" or $a[3+$e] eq "C" or $a[3+$e] eq "-") {
			$new_line .= "\t".$a[3+$e];
			$num_calls++;
			$num_calls_with_svdel++;
			if ($a[3+$e] ne $ref) {
				$num_snps++;
			}
		}
		else {
			if (defined($DELETIONS{$ECOTYPES[$e]}{$chr."#".$pos})) {
				$new_line .= "\tD";
				$num_calls_with_svdel++;
			}
			else {
				if (defined($UNSEQUENCED{$ECOTYPES[$e]}{$chr."#".$pos})) {
					$new_line .= "\tZ";
				}
				else {
					$new_line .= "\t".$a[3+$e];
				}
			}
		}
	}

	$new_line .= "\n"; # .$a[83]."\n"; # freq # Jun has not added this value anymore

	$NUM_CALLS{$num_calls}++;	
	$NUM_CALLS_WITH_SVDEL{$num_calls_with_svdel}++;

	if ($num_snps > 0) {
		$NUM_SNPS++;
		if ($num_snps > 1) {
			$NUM_SNPS_NONUNIQ++;
		}
	}
	
	print OUT $new_line;
}

# Print stats:
open STATS, ">stats_genome_matrix.$min_qual_value.txt";

print STATS "Number of calls (excl/incl SVdels):\n";
for (my $i = 0; $i <= 80; $i++) {
	print STATS $i, "\t", $NUM_CALLS{$i}, "\t", $NUM_CALLS_WITH_SVDEL{$i}, "\n";
}

print STATS "Number of SNPs: ", $NUM_SNPS, "\n";
print STATS "Number of SNPs (non-unique): ", $NUM_SNPS_NONUNIQ, "\n";



