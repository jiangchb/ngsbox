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
#  Module: Analysis::Eighties::THA_pipeline::run_THA_80.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 TargetRegionFile(absolute path required)\n";

my $THA_target_file = shift or die $usage;
open TF, $THA_target_file or die $usage;
my $target_num = 0;
while (<TF>) { $target_num++; }
close TF;

chdir("/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2");
my $base = "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2";
my $base_runfolder = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/EIGHTY";

my $shore_bin = "~/shore/startshore.sh";
my $live = 1;

my @ECOTYPES = ();

@ECOTYPES = ("Agu-1","Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE1","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE173","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169","ICE181","ICE21","ICE212","ICE213","ICE216","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE7","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2","Ped-0","Pra-6","Qui-0","Rue3-1","Sha","Star-8","TueSB30","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1");


print "Number of accession to be parsed:", @ECOTYPES+0, "\n";

for (my $i = 0; $i <@ECOTYPES; $i++) {
	my $ecotype = $ECOTYPES[$i];
	print "###########################################\n";
	print $ecotype, "\n";

	my $ecotype_folder = "/".$base."/".$ecotype."/AlignmentFolder";
	chdir($ecotype_folder);

	# Create library.txt file -- if not existent
	if (not -e "libraries.txt") {
		system("echo '1\tPE\t180.0\t400\t36' > libraries.txt");
	}	

	# Set map.list correctly
	my $maplist = "";
	my $tmpfolder = "";
	
	if (-e "map.list") {
		$maplist = "map.list";
        }
        else {
		if (not (-e "map.list.gz")) {
			die("No map.list(.gz) found for ".$ecotype."\n");
		}
        	$tmpfolder = "/ebio/abt6_projects2/nobackup/korbinian/EIGHTIES/".$ecotype;
                if (not (-e $tmpfolder)) {
			if ($live == 1) { mkdir($tmpfolder); }
			else { print "mkdir $tmpfolder\n"; }
                }
                $maplist = $tmpfolder."/map.list";
                if (not (-e $maplist)) {
			my $cmd_gzip = "gzip -cd map.list.gz > $maplist";
			if ($live == 1) { system($cmd_gzip); }
			else { print $cmd_gzip."\n";}
			
                }
        }

	# Merge left over
	if (not -e "left_over") {
		# collect run folder
		my $ecotype_runfolder_1 = "/".$base_runfolder."/".$ecotype."/run*";
		my $ecotype_runfolder_2 = "/".$base."/".$ecotype."/run*";
		my @runfolder = glob($ecotype_runfolder_1);
		my @runfolder2 = glob($ecotype_runfolder_2);

		for (my $j = 0; $j < @runfolder2; $j++) {
			push @runfolder, $runfolder2[$j];
		}

		# call shore mergeleftover
		my $cmd_mlo = "$shore_bin mergeleftover -p ".join(",",@runfolder)." -m $maplist -d left_over";

		if ($live == 1) { system($cmd_mlo); }
		else { print($cmd_mlo."\n"); }
		print "--------------------------------\n\n";
		print join(",",@runfolder), "\n";
		print "--------------------------------\n\n";	
	}

	# Call THA
	my $output_folder = "THA";
	if (not (-e $output_folder)) {
		my $cmd_shore = "$shore_bin THA -Z -L libraries.txt -i $maplist -l left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o $output_folder -v -y -M 36 -j 15,17,19,21 -T $THA_target_file";
		if ($live == 1) { system($cmd_shore); }
                else { print($cmd_shore."\n"); }

		for (my $num = 1; $num <= $target_num; $num++) {
			chdir($output_folder."/BuildingSite/Contigs/contig_".$num);
			my $cmd_amos = "AMOScmp contigs -D REF=../../Input/refseq_".$num.".fa";
			if ($live == 1) { system($cmd_amos); }
                	else { print($cmd_amos."\n"); }
			chdir("../../../..");
		}
	}


	if ($tmpfolder ne "") {
		my $cmd_rm = "nohup rm -rf $tmpfolder &";
		if ($live == 1) { system($cmd_rm); }
	        else { print($cmd_rm."\n"); }
	}

	

}


