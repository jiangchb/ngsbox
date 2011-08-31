#!/usr/bin/perl
use strict;
use warnings;

my $mapfile     = shift;
my $speciesfile = shift;


### Read species names
my $counter = 1;
my %species = ();
open IN, $speciesfile or die "Cannot open species file\n";
while ( <IN> ) {
	chomp;
	$species{$counter} = $_;
	$counter++;
}
close IN;


### Set species name in mapfile
open IN, $mapfile or die "Cannot open map file\n";
while ( <IN> ) {
	chomp;
	my @elem = split(/\t/, $_);

	print $species{$elem[0]} . "\t$elem[1]\t$elem[2]\t$elem[3]\t$elem[4]\t$elem[5]\t$elem[6]\n";
}
close IN;

exit(0);
