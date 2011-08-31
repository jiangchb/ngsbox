#! /usr/bin/perl
use strict;
my $file = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/Col-0-Geneva/AlignmentFolder/CNV_per_segment/shore_count.ALL_80.Bak-7.notcomplete.txt";
open FILE, $file or die "Cannot open file:".$file."\n";
open BACK, ">/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/Col-0-Geneva/AlignmentFolder/CNV_per_segment/shore_count.ALL_80.Bak-7.notcomplete.frac_back.txt";
open NBS, ">/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/Col-0-Geneva/AlignmentFolder/CNV_per_segment/shore_count.ALL_80.Bak-7.notcomplete.frac_nbs.txt";
while (my $line = <FILE>) {
	my @a = split " ", $line;
	for (my $i = 7; $i <= 246; $i+=3) {
		if ($i != 13) { # excl bak-7
		my $count = $a[$i];
		my $rpm = $a[$i+1];
		my $frac = $a[$i+2];
		if ($a[1] eq "NBS_LRR_active") {
			print NBS $frac, "\n"; 
		}
		print BACK $frac, "\n"; 
		}
	}	
}


