#! /usr/bin/perl

use strict;
use FindBin;

#chdir("/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/EIGHTY");
#my $base = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/EIGHTY/";

#my @ECOTYPES = ("Agu-1","Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE1","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE133","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169","ICE181","ICE21","ICE212","ICE213","ICE216","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE7","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2N8","Ped-0","Pra-6","Qui-0","Rue3.1","Sha","Star-8","Tue-SB30","Tue-SB300","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1");


chdir("/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/");
my $base = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/";

my @ECOTYPES = ("Col-0-Geneva");


for (my $i = 0; $i <@ECOTYPES; $i++) {
	my $ecotype = $ECOTYPES[$i];

	my $ecotype_folder = $base.$ecotype."/";
	my $inputfile = "$ecotype_folder/AlignmentFolder_DM2/DM2_coverage.txt";
	my $pdffile = "$ecotype_folder/AlignmentFolder_DM2/DM2_".$ecotype."_coverage.pdf";

	
	my $cut_cmd = "cut -f1,2,4,26 $ecotype_folder/AlignmentFolder_DM2/Analysis_DM2/ConsensusAnalysis/supplementary_data/consensus_summary.txt > $inputfile";
	print STDERR $cut_cmd, "\n";
	system($cut_cmd);

	my $r_cmd = "R --slave --vanilla --args $ecotype $inputfile $pdffile < ".$FindBin::Bin."/plot_coverage_DM2.R";
	print STDERR $r_cmd, "\n";
	system($r_cmd);

}


