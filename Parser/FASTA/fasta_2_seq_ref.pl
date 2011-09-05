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
#  Module: Parser::FASTA::fasta_2_seq_ref.pl
#  Purpose:
#  In:
#  Out:
#

# written by korbinian schneeberger


use Getopt::Long;

my %CMD;
my $assem;
my $file;
my $chr;
my $chr_count = 0;
my $pos = 1;

GetCom();

open FILE, $CMD{file} or die "Cannot open file:", $CMD{file},"\n";
open REF, ">".$CMD{file}.".seq_ref";
open MAX, ">".$CMD{file}.".seq_max";
open FAS, ">".$CMD{file}.".fa";

RUN: while (my $l = <FILE>) {
	chomp($l);
	if (substr($l, 0, 1) eq ">") {
		if (defined($chr)) {
			print MAX $assem, "\t", $chr_count, "\t", $chr, "\t", "\\N", "\t", $pos, "\n"; 
		}
		my @a = split " ", $l;
		$chr = substr($a[0], 1, length($a[0])-1);
		$chr_count++;
		$pos = 1;
		print FAS ">$chr_count\n";
	} else {
		if (defined($chr)) {
			for (my $i=0; $i<length($l); $i++) {
				print REF $chr_count, "\t", $pos, "\t", substr($l, $i, 1), "\n";
				$pos++;
			}
			print FAS "$l\n";
		}
	}
}
if (defined($chr)) {
	print MAX $assem, "\t", $chr_count, "\t", $chr, "\t", "\\N", "\t", $pos, "\n";
}

close FILE; close REF; close MAX; close FAS;
exit(0);

sub GetCom{

	my @usage = ("\nUsage: $0
                --file=file\treference genome
		--assem=name\tassembly name

	  	description:
		Makes a seq_ref and seq_max table from a fasta file
		\n\n");
 
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "file=s", "assem=s");

	die("Please specify file\n")	if not defined $CMD{file};
	die("Please specify assem\n")    if not defined $CMD{assem};

	$file = $CMD{file};
	$assem = $CMD{assem};

}
	      

