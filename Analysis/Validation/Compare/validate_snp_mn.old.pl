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
#  Module: Analysis::Validation::Compare::validate_snp_mn.old.pl
#  Purpose:
#  In:
#  Out:
#


my $mn_frag_trimming = 0;


my $usage = "$0 mn_fragments.txt mn_snps.txt shore_snps.txt offset\n";
my $frag = shift or die $usage;
my $mnsnps = shift or die $usage;
my $shoresnps = shift or die $usage;
my $offset = shift;

open FILE, $frag or die $usage;
my %POS = ();
while (my $line = <FILE>) {
	my @a = split " ", $line;
	for (my $i = $a[1]+$mn_frag_trimming; $i <= $a[2]-$mn_frag_trimming; $i++) {
		$POS{$a[0]."#".$i} = 1;
	}
}
close FILE;

open FILE, $mnsnps or die $usage;
my %MNSNPS = ();
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[0] <= 5) {
		if (defined($POS{$a[0]."#".$a[1]})) {
			$MNSNPS{$a[0]."#".$a[1]} = $a[2]."\t".$a[3];
		}
	}
}
close FILE;

open FILE, $shoresnps or die $usage;

my $tp = 0;
my $fp = 0;
my $fn = 0;

my %FALSE = ();

##################################################################################################################
##################################################################################################################

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[1] <= 5 && defined($POS{$a[1]."#".$a[2]})) {
		if (defined($MNSNPS{$a[1]."#".$a[2]})) {
			$tp++;
			delete $MNSNPS{$a[1]."#".$a[2]};
		}
		else {
			my $found = 0;
			for (my $i = 0; $i<$offset; $i++) {
				if (defined($MNSNPS{$a[1]."#".($a[2]-$i)}) or defined($MNSNPS{$a[1]."#".($a[2]-$i)})) {
					$found = 1;
				}
			}
			if ($found == 1) {
				$tp++;
				delete $MNSNPS{$a[1]."#".$a[2]};
			}
			else {
				$fp++;
				$FALSE{$a[1]."#".$a[2]} = "FP"."\t".$a[3]."\t".$a[4];
			}
		}
	}
}

foreach my $key (keys %MNSNPS) {
	$fn++;
	$FALSE{$key} = "FN\t".$MNSNPS{$key};
}


print STDERR "True Positives:", $tp, "\n";
print STDERR "False Positives:", $fp, "\n";
print STDERR "False Negatives:", $fn, "\n";

foreach my $key (sort keys %FALSE) {
	my @a = split "#", $key;
	print $a[0], "\t", $a[1], "\t", $FALSE{$key}, "\n";
}










