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
#  Module: Parser::SHORE::Variants::select_by_quality.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "\n$0 column homerAnno1 homerAnno2\n\n";
my $col   = shift or die $usage;
my $file1 = shift or die $usage;
my $file2 = shift or die $usage;

my %hash1 = ();
my %hash2 = ();


### Read merge column of file 1
open FILE1, $file1 or die $usage;
while(<FILE1>) {
	chomp;
	my @a = split(/\t/, $_);
	$hash1{$a[$col]} = 1;
}
close FILE1;


### Read merge column of file 2
open FILE2, $file2 or die $usage;
while(<FILE2>) {
	chomp;
	my @a = split(/\t/, $_);
	$hash2{$a[$col]} = 1;
}
close FILE2;


### Print new file 1
open FILE1, $file1 or die $usage;
open OUT1, ">$file1.merged" or die $usage;
while (<FILE1>) {

	chomp;
	my @a = split(/\t/, $_);

	if(exists $hash2{$a[$col]} ) {
		print OUT1 "$_\tYES\n";
	}
	else {
		print OUT1 "$_\tNO\n";
	}
}
close FILE1;
close OUT1;

### Print new file 2
open FILE2, $file2 or die $usage;
open OUT2, ">$file2.merged" or die $usage;
while (<FILE2>) {
	chomp;
	my @a = split(/\t/, $_);

	if(exists $hash1{$a[$col]} ) {
		print OUT2"$_\tYES\n";
	}
	else {
		print OUT2"$_\tNO\n";
	}
}
close FILE2;
close OUT2;

exit(0);

