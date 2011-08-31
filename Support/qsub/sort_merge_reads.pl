#! /usr/bin/perl
use strict;
use warnings;

my $usage = "\n$0 infolders\n\n";
my $ins  = shift or die $usage;


my @infolders = split(",", $ins);
my $files_r1 = "";
my $files_r2 = "";


foreach my $infolder (@infolders) {

	### Dir Level 1
	my @typefolders = glob($infolder . "/*");

	foreach my $typefolder (@typefolders) {

		### Dir Level 2
		my @typepath = split("/", $typefolder);
		my $typeleaf = $typepath[$#typepath];

		my @lengthfolders = glob($typefolder . "/*");

		foreach my $lengthfolder (@lengthfolders) {
		
			### Dir Level 3
			my @lengthpath = split("/", $lengthfolder);
			my $lengthleaf = $lengthpath[$#lengthpath];

			if($typeleaf eq "1") {
				if(-e "$lengthfolder/reads_0.fl") {
					$files_r1 .= "$lengthfolder/reads_0.fl,"
				}
				elsif(-e "$lengthfolder/reads_0.fl.gz") {
					$files_r1 .= "$lengthfolder/reads_0.fl.gz,"
				}
			}
			elsif($typeleaf eq "2") {
				if(-e "$lengthfolder/reads_0.fl") {
					$files_r2 .= "$lengthfolder/reads_0.fl,"
				}
				elsif(-e "$lengthfolder/reads_0.fl.gz") {
					$files_r2 .= "$lengthfolder/reads_0.fl.gz,"
				}
			}
		}	
	}
}

chop($files_r1);
chop($files_r2);
print "\n\nshore sort -o LongIndels/reads_1.fl -m -k 1i -i $files_r1\n\n";
print "\n\nshore sort -o LongIndels/reads_2.fl -m -k 1i -i $files_r2\n\n";

