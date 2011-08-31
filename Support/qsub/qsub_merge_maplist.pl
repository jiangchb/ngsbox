#! /usr/bin/perl
use strict;
use warnings;

my $usage = "\n$0 runfolder outfolder qsubname\n\n";
my $runfolder  = shift or die $usage;
my $qsub_name  = shift or die $usage;

open OUT, ">$qsub_name" or die $usage;


### Print qsub header and exports
print OUT "#!/bin/bash\n\n";
print OUT "#\$ -e $runfolder\n";
print OUT "#\$ -o $runfolder\n\n";
print OUT "source /users/GD/so/sossowski/.bashrc\n";
print OUT "export TMPDIR=/users/GD/projects/familiaMarruecos/tmp/\n\n";


my $counter = 1001;

my @lanefolders = glob($runfolder . "/*");

foreach my $lanefolder (@lanefolders) {

	my @typefolders = glob($lanefolder . "/*");

	foreach my $typefolder (@typefolders) {
		
		my @subpath = split("/", $typefolder);
		my $subleaf = $subpath[$#subpath];

		my @lengthfolders = glob($typefolder . "/*");

		foreach my $lengthfolder (@lengthfolders) {
			print OUT "mv $lengthfolder/map.list.gz $lengthfolder/map.first.gz\n\n";
			print OUT "mv $lengthfolder/map.list.gzx.gz $lengthfolder/map.first.gzx.gz\n\n";

			if($subleaf ne "single") {
				print OUT "shore sort -m -k4i -o $lengthfolder/map.list.gz $lengthfolder/map.first.gz $lengthfolder/map.left.gz\n\n";
			}
			else {
				print OUT "shore sort -m -k1i2i -o $lengthfolder/map.list.gz $lengthfolder/map.first.gz $lengthfolder/map.left.gz\n\n\n\n";
			}
		}
	}
}
