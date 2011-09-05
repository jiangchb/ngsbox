#! /usr/bin/perl
use strict;
use warnings;

###### 
# NGSbox - bioinformatics analysis tools for next generation sequencing data
#
# Copyright 2007-2011 Stephan Ossowski, Korbinian Schneeberger
# 
# NGSbox is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or any later version.
#
# NGSbox is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# Please find the GNU General Public License at <http://www.gnu.org/licenses/>.
#
#  -------------------------------------------------------------------------
#
#  Module: Analysis::Metagenomics::parse_species_and_reads.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "perl parse_species_and_reads.pl map.final\n";
my $file = shift or die $usage;

open FILE, $file;

open TMP, ">tmp.file";

while (<FILE>) {
	my @a = split "\t";
	my @b = split " ", $a[0];
	print TMP $a[3], "\t", $b[1], " ", $b[2], "\n"; 
}

close TMP;

system("sort tmp.file > tmp.file.sorted");

open TMP, "tmp.file.sorted";
open OUT, ">reads_by_species.txt";
my $curr = "";
my @s = ();
my %UNIQ_BY_SPECIES = ();

while (<TMP>) {
	chomp();
	my @a = split "\t";
	if ($curr != "" and $curr != $a[0]) {
		print OUT $curr, "\t";
		my %dup = ();
		foreach my $spe (sort @s) {
			if (not defined ($dup{$spe})) {
				print OUT $spe, ",";
				$UNIQ_BY_SPECIES{$spe}++;
			}
			$dup{$spe} = 1;
		}
		print OUT "\n";
		@s = ();	
	}
	push @s, $a[1];
	$curr = $a[0];
}

print OUT $curr, "\t";
foreach my $spe (sort @s) {
	print OUT $spe, "\,";
}
print OUT "\n";

close OUT;

system("cut -f2 reads_by_species.txt| sort | uniq -c | sort -n > reads_count_species.txt");
system("rm tmp.file tmp.file.sorted");

open OUT, ">tmp.file";

foreach my $spec (keys %UNIQ_BY_SPECIES) {
	print OUT $spec, "\t", $UNIQ_BY_SPECIES{$spec}, "\n";
}

close OUT;

system("sort -k3n tmp.file > reads_in_species.txt ");
system("rm tmp.file");

