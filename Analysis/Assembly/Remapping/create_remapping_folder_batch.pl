#! /usr/bin/perl
use strict;

my $usage = "$0 remappingfolder(will_be_created) runfolder runfolder ...\n";

my $remapping_folder = shift or die $usage;

my @reads_folder_names = ("1", "2", "single");


if (-e $remapping_folder) {
	print STDERR "Remapping folder exists already.\n";
	exit(1);
}

mkdir($remapping_folder) or die "Cannot create folder\n";

for (my $i = 0; $i < @ARGV; $i++) {

	my @a = split "/", $ARGV[$i];
	my $folder = $a[$#a];
	print "\n$folder\n";

	mkdir($remapping_folder."/".$folder);

	my @runfolders = glob($folder."/run*");

	# Parse lanes
	foreach my $runfolder (@runfolders) {

		my @b = split "/", $runfolder;
		my $run = $b[$#b];
		mkdir($remapping_folder."/".$folder."/".$run);
		print $run."\n";

		for (my $lane = 1; $lane <= 8; $lane++) {
			
			if (-e $ARGV[$i]."/".$run."/".$lane) {
				print $remapping_folder."/".$folder."/".$run."/".$lane."\n";
				mkdir($remapping_folder."/".$folder."/".$run."/".$lane);

				# Parse read folder
				for (my $read = 0; $read <= 2; $read++) {
	
					my $read_folder = $reads_folder_names[$read];
	
					if (-e $ARGV[$i]."/".$run."/".$lane."/".$read_folder) {
						mkdir($remapping_folder."/".$folder."/".$run."/".$lane."/".$read_folder);
		
						# Parse Length Folder
						my @length_folder = glob($ARGV[$i]."/".$run."/".$lane."/".$read_folder."/length_*");
	
						for (my $length = 0; $length < @length_folder; $length++) {
							my @b = split "/", $length_folder[$length];
							my $lfolder = $b[$#b];
	
							mkdir($remapping_folder."/".$folder."/".$run."/".$lane."/".$read_folder."/".$lfolder);
	
							### Copy reads_0.fl files and gzip
							if( -e $ARGV[$i]."/".$run."/".$lane."/".$read_folder."/".$lfolder."/reads_0.fl" ) {
								system("gzip -9 ".$ARGV[$i]."/".$run."/".$lane."/".$read_folder."/".$lfolder."/reads_0.fl");
							}
							system("cp ".$ARGV[$i]."/".$run."/".$lane."/".$read_folder."/".$lfolder."/reads_0.fl.gz ".$remapping_folder."/".$folder."/".$run."/".$lane."/".$read_folder."/".$lfolder."/.");
						}
					}
				}
			}
		}
	}
}

