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
#  Module: Simulation::Simulate_NGS::generate_meth_conv_genome.pl
#  Purpose:
#  In:
#  Out:
#



my $file = shift;
my $reffile = shift;

open FILE, $file;

my %CTPOS = ();
my %GAPOS = ();
while (<FILE>) {
	my @a = split " ";
	if ($a[3] > $a[2] and $a[3] > 5) {
#print $a[0], "#", $a[1], "\n";
		$CTPOS{$a[0]."#".$a[1]} = 1;
	}	
	if ($a[5] >= 5 and $a[5] > $a[4]) {
                $GAPOS{$a[0]."#".$a[1]} = 1;
        }
}

close FILE;

open FILE, $reffile;
my $chr = 0;
my $pos = 0;

open GAFILE, ">conv.ga.fa";
open CTFILE, ">conv.ct.fa";

my $count_ct = 0;
my $count_ga = 0;

while (my $line = <FILE>) {
	chomp($line);
	my @a = split " ", $line;
	if (substr($line, 0, 1) eq ">") {
		$chr = substr($a[0], 1, length($a[0]) -1);
		print STDERR $chr, "\n";
		print GAFILE $line, "\n";
		print CTFILE $line, "\n";
	}
	else {
		for (my $i = 0; $i < length($line); $i++) {
			$pos++;
			if (substr($line, $i, 1) eq "C" and defined($CTPOS{$chr."#".$pos})) {
				print CTFILE "T";
				print GAFILE substr($line, $i, 1);
				$count_ct++;
			}
			elsif (substr($line, $i, 1) eq "G" and defined($GAPOS{$chr."#".$pos})) {
				print GAFILE "A";
                                print CTFILE substr($line, $i, 1);
				$count_ga++;
                        }
			else {
				print GAFILE substr($line, $i, 1);
				print CTFILE substr($line, $i, 1);
			}
			if ($pos % 60 == 0) {
				print GAFILE "\n";
				print CTFILE "\n";
			}
		}
	}
}

print STDERR "Converted CT:", $count_ct, "\n";
print STDERR "Converted GA:", $count_ga, "\n";


