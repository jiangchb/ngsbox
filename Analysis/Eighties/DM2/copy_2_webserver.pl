#! /usr/bin/perl

use strict;
use FindBin;

chdir("/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2");

my @ECOTYPES = ("Agu-1","Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE1","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE133","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169","ICE181","ICE21","ICE212","ICE213","ICE216","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE7","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2N8","Ped-0","Pra-6","Qui-0","Rue3.1","Sha","Star-8","Tue-SB30","Tue-SB300","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1");


for (my $i = 0; $i <@ECOTYPES; $i++) {
	my $ecotype = $ECOTYPES[$i];
	print STDERR $ecotype, "\n";
	my $cp_cmd = "cp $ecotype/AlignmentFolder/map.list.gz tmp/.";
	system($cp_cmd);
	my $mv_cmd = "mv tmp/map.list.gz tmp/$ecotype.map.list.gz";
        system($mv_cmd);

}


