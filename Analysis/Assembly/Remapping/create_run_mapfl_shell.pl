#! /usr/bin/perl
use strict;

my $usage = "$0 project_folder\n";

my $project_folder = shift or die $usage;




my @strain_folders = glob($project_folder."/*");

# Parse strain folders
my @strain_folders = glob($project_folder."/*");
foreach my $strain_folder (sort @strain_folders) {

	# Parse run folders
	my @run_folders = glob($strain_folder."/run_*");
	foreach my $run_folder (sort @run_folders) {

		# ~/shore_mapfl/shore mapflowcell -i /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked_All_Mapper/TAIR8.v1.Masked.fa.shore -f . -u -n 7% -g 1 -D 50000 -S -L 14 -R 1 -c 16 -b 200000 -v blat -p
		print "~/shore_mapfl/shore mapflowcell -i /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked_All_Mapper/TAIR8.v1.Masked.fa.shore -f $run_folder -u -n 7% -g 1 -D 1000 -S -L 14 -R 1 -c 16 -b 200000 -v blat -p\n";
	}
}

