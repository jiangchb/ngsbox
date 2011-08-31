#! /usr/bin/perl

use strict;

my $usage = "$0 batch\n";
my $batch = shift or die $usage;

chdir("/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2");
my $base = "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2";

my $shore_bin = "shore_rur";

my @ECOTYPES = ();

@ECOTYPES = ("Agu-1","Bak-2","Bak-7") if $batch == 1;
@ECOTYPES = ("Cdm-0","Del-10","Dog-4") if $batch == 2;
@ECOTYPES = ("Don-0","Ey15-2","Fei-0") if $batch == 3;
@ECOTYPES = ("HKT2.4","ICE1","ICE102","ICE104","ICE106") if $batch == 4;
@ECOTYPES = ("ICE107","ICE111","ICE112","ICE119","ICE120") if $batch == 5;
@ECOTYPES = ("ICE127","ICE130","ICE173","ICE134","ICE138") if $batch == 6;
@ECOTYPES = ("ICE150","ICE152","ICE153","ICE163","ICE169") if $batch == 7;
@ECOTYPES = ("ICE181","ICE21","ICE212","ICE213","ICE216") if $batch == 8;
@ECOTYPES = ("ICE226","ICE228","ICE29","ICE33","ICE36") if $batch == 9;
@ECOTYPES = ("ICE49","ICE50","ICE60","ICE61","ICE63") if $batch == 10;
@ECOTYPES = ("ICE7","ICE70","ICE71","ICE72","ICE73") if $batch == 11;
@ECOTYPES = ("ICE75","ICE79","ICE91","ICE92","ICE93") if $batch == 12;
@ECOTYPES = ("ICE97","ICE98","Istisu-1","Kastel-1","Koch-1") if $batch == 13;
@ECOTYPES = ("Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6") if $batch == 14;
@ECOTYPES = ("Nie1-2","Ped-0","Pra-6") if $batch == 15;
@ECOTYPES = ("Qui-0","Rue3-1") if $batch == 16;
@ECOTYPES = ("Sha","Star-8","TueSB30") if $batch == 17;
@ECOTYPES = ("Tuescha9","TueV12") if $batch == 18;
@ECOTYPES = ("TueWa1-2","Vash-1","Vie-0") if $batch == 19;
@ECOTYPES = ("WalhaesB4","Xan-1","Yeg-1") if $batch == 20;

# NICHT FUNKTIONIERT:

my @ECOTYPES = ("Bak-7") if $batch == 21;


print STDERR "Number of accession to be parsed:", @ECOTYPES+0, "\n";

for (my $i = 0; $i <@ECOTYPES; $i++) {
	my $ecotype = $ECOTYPES[$i];

	my $ecotype_folder = "/".$base."/".$ecotype."/";
	my $alignment_folder = $ecotype_folder."/AlignmentFolder";

	# Create new runfolders

	chdir($alignment_folder);
	print $alignment_folder, "\n";

	my $output_folder = "SV_lib1_version4_final";
	my $maplist = "";

	if (not (-e $output_folder)) {

		my $tmpfolder = "";

		if (-e "map.list") {	
			$maplist = "map.list";
		}
		else {
			$tmpfolder = "/ebio/abt6_projects2/nobackup/korbinian/EIGHTIES/".$ecotype;
			if (not (-e $tmpfolder)) {
				mkdir($tmpfolder);
			}

			$maplist = $tmpfolder."/map.list";
			if (not (-e $maplist)) {
				system("gzip -cd map.list.gz > $maplist");
			}
		}

		# Call shore structure all read pairs:
		system("/ebio/abt6/korbinian/shore/$shore_bin struct -E ids.false_mapping.ignore.lib1.txt -e 6,7 -N $ecotype -s 200 -u 0 -l 1 -T PE -i $maplist -o $output_folder -t /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore.trans -C /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.centromeres.txt -b Analysis_q01/ConsensusAnalysis/supplementary_data/orphan_distribution.txt -c Analysis_q01/ConsensusAnalysis/supplementary_data/unseq_core.txt -h Analysis_q01/ConsensusAnalysis/supplementary_data/unseq_cn.txt -O Analysis_q01/ConsensusAnalysis/supplementary_data/oversampled.txt");

		if ($tmpfolder ne "") {
			system("nohup rm -rf $tmpfolder &");
		}

	}

}


