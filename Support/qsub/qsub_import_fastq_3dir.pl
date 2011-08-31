#! /usr/bin/perl
use strict;
use warnings;

my $usage = "\n$0 infolder outfolder qsubName\n\n";
my $infolder  = shift or die $usage;
my $outfolder = shift or die $usage;
my $qsub_name = shift or die $usage;

open OUT, ">$qsub_name" or die $usage;


### Print qsub header and exports
print OUT "#!/bin/bash\n\n";
print OUT "#\$ -e $outfolder\n";
print OUT "#\$ -o $outfolder\n\n";
print OUT "export SHORE=/users/GD/tools/shore/\n";
print OUT "export TMPDIR=/users/GD/projects/familiaMarruecos/tmp/\n\n";

my $counter = 1007;

### Dir Level 1
my @samplefolders = glob($infolder . "/*");

foreach my $samplefolder (@samplefolders) {

	### Dir Level 2
	my @samplepath = split("/", $samplefolder);
	my $sampleleaf = $samplepath[$#samplepath];

	my @FCfolders = glob($samplefolder . "/*");

	#if(! -e "$outfolder/$sampleleaf") {
		#mkdir("$outfolder/$sampleleaf");
	#}

	foreach my $FCfolder (@FCfolders) {
		
		### Dir Level 3
		my @FCpath = split("/", $FCfolder);
		my $FCleaf = $FCpath[$#FCpath];

		if($FCleaf =~ /^s_/) {
			next;
		}

		my @files = glob($FCfolder . "/*");

		my @read1 = ();
		my @read2 = ();

		foreach my $file (@files) {
			my @filepath = split("/", $file);
			my $fileleaf = $filepath[$#filepath];
	
			my @junk = split("_", $fileleaf);
		
			if($junk[2] == 1) {
				push @read1, $file;
			}
			else {
				push @read2, $file;
			}
		}
	
		print OUT "/users/GD/tools/shore/shore import -v Fastq -e Shore -a genomic -i $counter -o $outfolder/$sampleleaf/$FCleaf -n 4 -g -c -V 10 -x ";
		print OUT join ",", @read1;
		print OUT " -y ";
		print OUT join ",", @read2;
		print OUT "\n\n";
	
		$counter++;
	}
}
