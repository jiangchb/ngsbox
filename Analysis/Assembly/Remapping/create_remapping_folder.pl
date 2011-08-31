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

	mkdir($remapping_folder."/".$folder);

	# Parse lanes
	for (my $lane = 1; $lane <= 8; $lane++) {
		if (-e $ARGV[$i]."/".$lane) {
			mkdir($remapping_folder."/".$folder."/".$lane);

			# Parse read folder
			for (my $read = 0; $read <= 2; $read++) {

				my $read_folder = $reads_folder_names[$read];

				if (-e $ARGV[$i]."/".$lane."/".$read_folder) {
					mkdir($remapping_folder."/".$folder."/".$lane."/".$read_folder);
	
					# Parse Length Folder
					my @length_folder = glob($ARGV[$i]."/".$lane."/".$read_folder."/length_*");

					for (my $length = 0; $length < @length_folder; $length++) {
						my @b = split "/", $length_folder[$length];
						my $lfolder = $b[$#b];

						mkdir($remapping_folder."/".$folder."/".$lane."/".$read_folder."/".$lfolder);

						#system("cp -l ".$ARGV[$i]."/".$lane."/".$read_folder."/".$lfolder."/left_over.blt ".$remapping_folder."/".$folder."/".$lane."/".$read_folder."/".$lfolder."/.");
						#system("cp -l ".$ARGV[$i]."/".$lane."/".$read_folder."/".$lfolder."/reads_0.fl ".$remapping_folder."/".$folder."/".$lane."/".$read_folder."/".$lfolder."/.");
						
						### Copy reads_0.fl files and gzip
						if( -e $ARGV[$i]. "/". $lane ."/". $read_folder ."/". $lfolder. "/reads_0.fl" ) {
							system("cp ".$ARGV[$i]."/".$lane."/".$read_folder."/".$lfolder."/reads_0.fl ".$remapping_folder."/".$folder."/".$lane."/".$read_folder."/".$lfolder."/.");
							system("gzip -9 ".$ARGV[$i]."/".$lane."/".$read_folder."/".$lfolder."/reads_0.fl ".$remapping_folder."/".$folder."/".$lane."/".$read_folder."/".$lfolder."/reads_0.fl");
						}
						elsif( -e $ARGV[$i]. "/". $lane ."/". $read_folder ."/". $lfolder. "/reads_0.fl.gz" ) {
							system("cp ".$ARGV[$i]."/".$lane."/".$read_folder."/".$lfolder."/reads_0.fl.gz ".$remapping_folder."/".$folder."/".$lane."/".$read_folder."/".$lfolder."/.");
						}
					}
				}
			}
		}
	}
}


