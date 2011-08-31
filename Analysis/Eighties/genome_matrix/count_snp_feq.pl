#! /usr/bin/perl
##################################################################################
# count snp includes deletion fre of each position in genome_matrix.base_call file
##################################################################################


use strict;

my $usage = "$0 genomeMatrix.base_call\n";

my $file = shift or die $usage;

open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;

	my $freq = 0;

	for (my $i = 3; $i < @a; $i+=1) {
			if ($a[$i] ne $a[2] and $a[$i] ne 'N' ) {
				$freq++;

		}

	}

	chomp($line);
        print $freq, "\n";
#	print $line, "\t", $freq, "\n";

}


