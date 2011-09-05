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
#  Module: Parser::FASTA::sort_fasta_header_human.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 fastafile\n";
my $file = shift or die $usage;
open FILE, $file or die $usage;

my %entries = ();
my $id = "";

while (my $line = <FILE>) {

	if (substr($line, 0, 1) eq ">") {
		my @e = split(" ", $line);
		$id = substr($e[0], 1);

		# Human specific modifications
		if($id eq "X")     { $id = 23; }
		elsif($id eq "Y")  { $id = 24; }
		elsif($id eq "MT") { $id = 25; }

		$entries{$id} = $line;
	}
	else {
		$entries{$id} .= $line;
	}
}

foreach my $header ( sort {$a<=>$b} keys %entries ) {
	print $entries{$header};
}
	

