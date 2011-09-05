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
#  Module: Analysis::Validation::Compare::get_ecotype_sequence.pl
#  Purpose:
#  In:
#  Out:
#


use lib "$ENV{PGSP}/Prediction/Validation/";
use Genome;

my $usage = "$0 chrfasta snpfile delfile insfile chr begin end\n";

my $reffile = shift or die $usage;
my $snpfile = shift or die $usage;
my $deletionfile = shift or die $usage;
my $insertionfile = shift or die $usage;
my $chr = shift or die $usage;
my $begin = shift or die $usage;
my $end = shift or die $usage;

my $genome = new Genome();
$genome->calc_chr_seq($reffile, $snpfile, $deletionfile, $insertionfile);

my $seq = $genome->get_sequence($chr, $begin, $end);

print ">", $chr, " ", $begin, " ", $end, "\n";
for (my $i = 1; $i<=length($seq); $i++) {
	print substr($seq, $i-1, 1);
	print "\n" if $i%60 == 0;
}




