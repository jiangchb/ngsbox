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
#  Module: Analysis::Eighties::CNV::call_cnv.pl
#  Purpose:
#  In:
#  Out:
#




chdir("/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2");
my $base = "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2";

my $shore_bin = "shore_ohm";
my $anno_file = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/Col-0-Geneva/AlignmentFolder/CNV_per_segment/TAIR8_transposable_element_gene.annotation.txt";

my @ECOTYPES = ("Agu-1","Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE1","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE173","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169","ICE181","ICE21","ICE212","ICE213","ICE216","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE7","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2","Ped-0","Pra-6","Qui-0","Rue3-1","Sha","Star-8","Tue-SB30","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1");
#my @ECOTYPES = ("Del-10","ICE138","ICE163","ICE228","ICE50","ICE73","TueSB30");
#my @ECOTYPES = ("Lag2.2");

print STDERR "Number of accession to be parsed:", @ECOTYPES+0, "\n";

for (my $i = 0; $i <@ECOTYPES; $i++) {
	my $ecotype = $ECOTYPES[$i];

	my $ecotype_folder = "/".$base."/".$ecotype."/";
	my $alignment_folder = $ecotype_folder."/AlignmentFolder";

	my $maplist = $alignment_folder."/map.list";
	my $tmpfolder = "";

	print $maplist, "\n";

	if (not (-e $maplist)) {	
		print "TMP\n";
		$tmpfolder = "/ebio/abt6_projects2/nobackup/korbinian/EIGHTIES_CNV/".$ecotype;
		if (not (-e $tmpfolder)) {
			mkdir($tmpfolder);
		}

		$maplist = $tmpfolder."/map.list";
		if (not (-e $maplist)) {
			system("gzip -cd $alignment_folder/map.list.gz > $maplist");
		}
	}

	my $output = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/Col-0-Geneva/AlignmentFolder/CNV_per_segment/shore_count.TE.".$ecotype.".txt";

	# Call shore structure all read pairs:
	
	my $cmd = "/ebio/abt6/korbinian/shore/$shore_bin count -m /ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/Col-0-Geneva/AlignmentFolder/map.list,$maplist -o $output -f $anno_file -W const";
	#print $cmd, "\n";
	system($cmd);
	
	if ($tmpfolder ne "") {
		system("rm -rf $tmpfolder");
	}

}


