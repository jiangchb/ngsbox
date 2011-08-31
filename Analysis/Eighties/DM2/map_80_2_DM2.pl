#! /usr/bin/perl

use strict;

chdir("/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/EIGHTY");
my $base = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/EIGHTY/";

my @ECOTYPES = ("Agu-1","Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE1","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE133","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169","ICE181","ICE21","ICE212","ICE213","ICE216","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE7","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2N8","Ped-0","Pra-6","Qui-0","Rue3.1","Sha","Star-8","Tue-SB30","Tue-SB300","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1");


#chdir("/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/");
#my $base = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/";


#my @ECOTYPES = ("Col-0-Geneva");



for (my $i = 0; $i <@ECOTYPES; $i++) {
	my $ecotype = $ECOTYPES[$i];

	my $ecotype_folder = $base.$ecotype."/";

	# Create new runfolders

	chdir($ecotype_folder);
	print $ecotype, "\n";

	if (not (-e "AlignmentFolder_DM2")) {

		my @runfolders = glob("run*");
		my @mappedfolders = ();
	
		foreach my $r (@runfolders) {
			my $runfolder = $ecotype_folder.$r;
			my $runfolder_new = $runfolder."_DM2/";
			$runfolder .= "/";
	
			# Create folder and new reads files
	
			mkdir($runfolder_new);
			chdir($runfolder);
	
	
			my @lanefolders = glob("*");
	
			foreach my $l (@lanefolders) {
				if ($l eq "1" or $l eq "2" or $l eq "3" or $l eq "4" or $l eq "5" or $l eq "6" or $l eq "7" or $l eq "8") {
					my $lanefolder = $runfolder.$l."/";
					my $lanefolder_new = $runfolder_new.$l."/";
					mkdir($lanefolder_new);
	
					chdir($lanefolder);
					my @readfolders = glob("*");
	
					foreach my $r (@readfolders) {
						if ($r eq "1" or $r eq "2" or $r eq "single") {
							my $readfolder = $lanefolder.$r."/";
							my $readfolder_new = $lanefolder_new.$r."/";
							mkdir($readfolder_new);
	
							chdir($readfolder);
							my @lengthfolders = glob("length_*");
	
							foreach my $e (@lengthfolders) {
								my $lengthfolder = $readfolder.$e."/";
								my $lengthfolder_new = $readfolder_new.$e."/";
								mkdir($lengthfolder_new);
	
									system("cp $lengthfolder/reads_0.fl $lengthfolder_new/.");
							}
	
						}							
					}
				} 
			}
	
			print STDERR $runfolder_new, "\n";
	
			# Call shore map 
			system("/ebio/abt6/korbinian/shore/shore_090731_dm2 mapflowcell -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/DM2/DM2.fa.shore -o $runfolder_new -n 10% -g 7% -c 8 -p");
	
			push @mappedfolders, $runfolder_new;
	
			chdir($runfolder);
	
			
			foreach my $l (@lanefolders) {
				if ($l eq "1" or $l eq "2" or $l eq "3" or $l eq "4" or $l eq "5" or $l eq "6" or $l eq "7" or $l eq "8") {
					my $lanefolder_new = $runfolder_new.$l."/";
					my $shore_corr_cmd = "/ebio/abt6/korbinian/shore/shore_090731_dm2 correct4pe -l $lanefolder_new ";
					if ($r !~ /MP/) {
						$shore_corr_cmd .= " -x 200 ";
					}	
					else {
						$shore_corr_cmd .= " -x 4000 -m ";
					}
					print STDERR $shore_corr_cmd, "\n";
					system($shore_corr_cmd);
				}
			}
	
	
		}
	
		chdir($ecotype_folder);
	
		my $shore_merge_cmd = "/ebio/abt6/korbinian/shore/shore_090731_dm2 merge -d AlignmentFolder_DM2 -p ";
		my $rm_cmd = "rm -rf ";
	
		foreach my $mf (@mappedfolders) {
			$shore_merge_cmd .= $mf.",";
			$rm_cmd .= $mf." ";
		}
	
		system($shore_merge_cmd);
		system($rm_cmd);
	
		my $shore_cons_cmd ="/ebio/abt6/korbinian/shore/shore_090731_dm2 consensus -n $ecotype -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/DM2/DM2.fa.shore -o AlignmentFolder_DM2/Analysis_DM2 -i AlignmentFolder_DM2/map.list -a ~/shore/Analysis/scoring_matrices/scoring_matrix_hom.txt -r -v -b 0.51";
	
		system($shore_cons_cmd);

	}

}


