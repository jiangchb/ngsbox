#! /usr/bin/perl
###########################################################################
# parse positions having at least one SNP from genome_matrix.base_call file 
###########################################################################


use strict;

my $usage = "$0 genomeMatrix_basse_call file\n";

my $file = shift or die $usage;

open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;


	for (my $i = 3; $i < @a; $i+=1) {
		if ($a[$i] ne $a[2] and $a[$i] ne "N") {
	

	chomp($line);
	print $line, "\n";
        last;
        } 
          
        }
}


