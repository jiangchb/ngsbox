#! /usr/bin/perl
use strict;

my $usage = "$0 project_folder\n";

my $project_folder = shift or die $usage;



my $count = 0;
my @strain_folders = glob($project_folder."/*");

# Parse strain folders
my @strain_folders = glob($project_folder."/*");
foreach my $strain_folder (sort @strain_folders) {


	# Parse run folders
	my @run_folders = glob($strain_folder."/run_*");
	foreach my $run_folder (sort @run_folders) {
		
		my @lane_folders = glob($run_folder."/*");
		foreach my $lane_folder (sort @lane_folders) {

			my @a = split("/", $lane_folder);
			my $lane = $a[$#a];

			if( $lane =~ /[1-9]/ && $lane !~ /_/) {
				open OUT, ">>run_corr4pe_$count.sh" or die;
				print OUT "~/shore_mapfl/shore correct4pe -l $lane_folder -x 200 -e 1 -p\n";
				close OUT;
			}
		}
	}
	$count++;
}

