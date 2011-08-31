#! /usr/bin/perl

use strict;

chdir("/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/EIGHTY");
my $base = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/EIGHTY";
my $base2 = "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2";

my @ECOTYPES = ("Agu-1","Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE134","ICE138","ICE150","ICE152","ICE153","ICE169","ICE173","ICE181","ICE1","ICE212","ICE213","ICE216","ICE21","ICE226","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE70","ICE71","ICE72","ICE75","ICE79","ICE7","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2","Ped-0","Pra-6","Qui-0","Rue3-1","Sha","Star-8","TueSB30","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1", "ICE163","Lag2.2","ICE228","ICE73");
#my @ECOTYPES = ("Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE134","ICE138","ICE150","ICE152","ICE153","ICE169","ICE173","ICE181","ICE1","ICE212","ICE213","ICE216","ICE21","ICE226","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE70","ICE71","ICE72","ICE75","ICE79","ICE7","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2","Ped-0","Pra-6","Qui-0","Rue3-1","Sha","Star-8","TueSB30","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1", "ICE163","Lag2.2","ICE228","ICE73");

#my @ECOTYPES = ("Bak-7","Fei-0","ICE173","Ped-0","Qui-0","Rue3-1","Tue-SB30");
#my @ECOTYPES = ("ICE181","Ped-0","Dog-4");


print STDERR "Number of accession to be parsed:", @ECOTYPES+0, "\n";

for (my $i = 0; $i <@ECOTYPES; $i++) {
	my $ecotype = $ECOTYPES[$i];
	
	print STDERR $ecotype, "\n";

	my $ecotype_folder = "/".$base."/".$ecotype."/";
	my $alignment_folder = $ecotype_folder."/LeftOvers";

	my $maplist_file = "/".$base2."/".$ecotype."/AlignmentFolder/map.list";
	if (! (-e $maplist_file)) {
		mkdir("/ebio/abt6_projects/nobackup/korbinian/EIGHTYLO/".$ecotype);
		system("gzip -cd /$base2/$ecotype/AlignmentFolder/map.list.gz > /ebio/abt6_projects/nobackup/korbinian/EIGHTYLO/$ecotype/map.list");
		$maplist_file =  "/ebio/abt6_projects/nobackup/korbinian/EIGHTYLO/$ecotype/map.list";		
	}

	chdir($ecotype_folder);

	my @runs = glob("run_*");
	my $run_string = join(",", @runs);

	# Call shore structure all read pairs:
	system("/ebio/abt6/korbinian/shore_git/shore/startshore.sh mergeleftover -p $run_string -d $alignment_folder -m $maplist_file");
	
	#print("/ebio/abt6/korbinian/shore/shore_final mergeleftover -p $run_string -d $alignment_folder\n");


	####################################################################
	
	#system("mv $ecotype_folder/LeftOvers/left_over/left_over_missing.fl  $ecotype_folder/$ecotype.left_over_missing.fl");
	#print("mv $ecotype_folder/LeftOvers/left_over/left_over_missing.fl  $ecotype_folder/$ecotype.left_over_missing.fl\n");
	#system("rm -rf $ecotype_folder/LeftOvers");
	#print("rm -rf $ecotype_folder/LeftOvers\n");

	if (-e "/ebio/abt6_projects/nobackup/korbinian/EIGHTYLO/$ecotype/map.list") {
		unlink("/ebio/abt6_projects/nobackup/korbinian/EIGHTYLO/$ecotype/map.list");
	}

}


