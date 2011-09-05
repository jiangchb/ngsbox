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
#  Module: Analysis::SV::Remapping::remap_discordant_singleaccession.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 batch\n";
my $batch_num = shift or die $usage;

my $live = 1;

#my $base = "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/";
my $base = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/";
#my $base = "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2";

chdir($base);
my @ECOTYPES= ();
@ECOTYPES = ("Agu-1") if $batch_num == 1;
@ECOTYPES = ("Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE1","ICE102","ICE104","ICE106") if $batch_num == 2;
@ECOTYPES = ("ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE133","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169") if $batch_num == 3;
@ECOTYPES = ("ICE181","ICE21","ICE212","ICE213","ICE216","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63") if $batch_num == 4;
@ECOTYPES = ("ICE7","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1") if $batch_num == 5;
@ECOTYPES = ("Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2","Ped-0","Pra-6","Qui-0","Rue3-1","Sha","Star-8","Tuescha9","TueV12") if $batch_num == 6;
@ECOTYPES = ("TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1") if $batch_num == 7;
@ECOTYPES = ("ICE36","ICE49","ICE50","ICE60","ICE61","ICE63") if $batch_num == 8;
@ECOTYPES = ("TueSB30") if $batch_num == 9;

# Wiederholen:

@ECOTYPES = ("Bak-7") if $batch_num == 10;

my $lib = 1;
my $expected = 200;
my $mate_pair = ""; # "-m"
#my $insert_dist = "insert_dist.lib$lib.txt";
my $insert_dist = "insert_dist.txt";


my $folder = "run_unhappy_lib".$lib;
my $shore = "shore_falseid";
my $cores = 3;

my $gunzip = 1;
my $parse = 1;
my $align = 1;
my $correct = 1;
my $rename = 1;
my $del = 1;

for (my $i = 0; $i <@ECOTYPES; $i++) {
	my $ecotype = $ECOTYPES[$i];
	print $ecotype, "\n";

	my $alignment_folder = $base."/".$ecotype."/AlignmentFolder/";

	chdir($alignment_folder) if $live == 1;

	my $tmpfolder = "";
	my $maplist = "map.list";

	# Gunzip if necessary
	if (not (-e "map.list") and $gunzip == 1 and $live == 1) {
		$tmpfolder = "/ebio/abt6_projects2/nobackup/korbinian/EIGHTIES_FALSE_MAPPING/".$ecotype;
                if (not (-e $tmpfolder)) {
	                mkdir($tmpfolder);
                }

                $maplist = $tmpfolder."/map.list";
                if (not (-e $maplist)) {
        	        system("gzip -cd map.list.gz > $maplist");
                }
		#my $gz_cmd = "gunzip map.list.gz";
		#print $gz_cmd, "\n";
		#system("gunzip map.list.gz") if $live == 1;
	}

	# Parse putative false mappings
	if ($parse == 1) {
		my $parse_cmd = "perl ".$ENV{PGSP}."/Analysis/SV/Remapping/parse_unhappy.pl $maplist $lib $folder";
		print $parse_cmd, "\n";
		system($parse_cmd) if $live == 1;
	}

	# Align
	if ($align == 1) {
		my $map_cmd = "/ebio/abt6/korbinian/shore/$shore mapflowcell -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o $folder -n 10% -g 7% -c $cores -a -p";
		print $map_cmd, "\n";
		system($map_cmd) if $live == 1;
	}
	
	# Correct 4 pairing
	if ($correct == 1) {
		my $correct4pe_cmd = "/ebio/abt6/korbinian/shore/$shore correct4pe -l $folder/1 -x $expected -e $lib -i $insert_dist $mate_pair";
		print $correct4pe_cmd, "\n";
		system($correct4pe_cmd) if $live == 1;
	}

	# parse happy 
	my $parse_happy_cmd = "perl ".$ENV{PGSP}."/Analysis/SV/Remapping/parse_happy.pl $folder/1/1/map.list.1 > ids.false_mapping.ignore.lib".$lib.".txt";
	print $parse_happy_cmd, "\n";
	system($parse_happy_cmd) if $live == 1;

	# del tmp folder
	if ($live == 1 and $tmpfolder ne "" and $del == 1) {
		my $del_cmd = "rm -rf $tmpfolder";
		print $del_cmd, "\n";
		system($del_cmd);
	}
	
	if ($live == 1 and $del == 1) {
		my $del_cmd2 = "rm -rf $folder";
                print $del_cmd2, "\n";
                system($del_cmd2);
	}

}


