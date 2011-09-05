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
#  Module: Analysis::CNV::anno2shore.pl
#  Purpose:
#  In:
#  Out:
#

my $usage = "$0 NBSLRR TAIR8_genes_trasposons shorecout\n";

my $nbslrr = shift or die $usage;
my $tair8 = shift or die $usage;
my $file = shift or die $usage;

#########################################################################################
my %NBS = ();
open FILE, $nbslrr or die "Cannot open file\n";
while (my $line = <FILE>) {
	chomp($line);
	$NBS{$line} = 1;
}
close FILE;

#########################################################################################
my %POS2ANNO = ();
open FILE, $tair8 or die "Cannot open file\n";
while (my $line = <FILE>) {
	if (substr($line, 0, 1) ne "#") {
		my @a = split " ", $line;
		$POS2ANNO{$a[1]."#".$a[2]."#".$a[3]} = $a[4];
	}
}
close FILE;

#########################################################################################
open FILE, $file or die "Cannot open file\n";

while (my $line = <FILE>) {
	if (not(substr($line, 0, 1) eq "#") and length($line) > 5)  {
		my @a = split " ", $line;

		chomp($line);
		print $line;

		if (defined($POS2ANNO{$a[1]."#".$a[2]."#".$a[3]})) {
			print "\t", $POS2ANNO{$a[1]."#".$a[2]."#".$a[3]};
			if (defined($NBS{$POS2ANNO{$a[1]."#".$a[2]."#".$a[3]}})) {
				print "\tNBS_LRR";
			}
			else {
				print "\tNULL";
			}
		}
		else {
			print "\tNULL\tNULL";
		}

		print "\n";
	}
}	

close FILE;



