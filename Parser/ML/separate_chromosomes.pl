#! /usr/bin/perl
use warnings;
use strict;

my $usage = "$0 maplist outfolder\n";

my $file      = shift or die $usage;
my $outfolder = shift or die $usage;

open FILE, $file or die "cannnot open $file\n";

my $chromosome = -1;

while(<FILE>) {
	my $line = $_;
	my @entries = split(/\t/, $line);
	if($entries[0] != $chromosome) {
		if($chromosome != -1) { close OUT; }
		$chromosome = $entries[0];
		open(OUT, ">$outfolder/map_chr$chromosome.list") or die "cannot open $outfolder/map_chr$chromosome.list\n";
	}

	print OUT $line;
}

exit(0);
