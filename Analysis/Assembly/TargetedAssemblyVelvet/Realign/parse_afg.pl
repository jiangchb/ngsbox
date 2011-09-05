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
#  Module: Analysis::Assembly::TargetedAssemblyVelvet::Realign::parse_afg.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "perl parse_agf.pl afg_file input_fasta\n";
open FILE, $ARGV[0] or die $usage;
open FASTA, $ARGV[1] or die $usage;

my $count = 1;
my %READ_ID = ();
while (my $line = <FASTA>) {
	chomp($line);
	if (substr($line, 0, 1) eq ">") {
		$READ_ID{$count} = substr($line, 1, length($line)-1);
		$count++;
	}
}

print STDERR "Parsed fasta\n";

my $id;
my $first_flag = 1;
CONTIGS: while (my $line = <FILE>) {
	chomp($line);
	if ($line eq "{CTG") {
		$id = <FILE>;
		chomp($id);
		$id =~ s/iid://g;
		if ($first_flag == 0) {
			print "\n";
		}
		else {
			$first_flag = 0;
		}
		print $id;
		<FILE>; <FILE>;
		my $seq = "";
		SEQ: while (my $l = <FILE>) {
			chomp($l);
			if (substr($l, 0, 1) eq ".") {
				print "\t", length($seq), "\t", $seq;
				last SEQ;	
			}
			$seq .= $l;	
		}
	}
	elsif ($line eq "{TLE") {
		$tile = <FILE>;
		chomp($tile);
		$tile =~ s/src://g;
		print "\t", $READ_ID{$tile};
	}
	elsif ($line eq "{RED") {
		last CONTIGS;
	}
}



