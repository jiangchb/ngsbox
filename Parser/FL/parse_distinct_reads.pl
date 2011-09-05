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
#  Module: Parser::FL::parse_distinct_reads.pl
#  Purpose:
#  In:
#  Out:
#

#written by ks

use Getopt::Long;

my %CMD;

GetCom();

my %IDS = ();
while (<FILE2>) {
	my @a = split " ";
	$IDS{$a[0]} = 1;
}


while (<FILE1>) {
        my @a = split " ";
	if (not defined $IDS{$a[0]}) {
		print $_;
	}
}

sub GetCom{

  my @usage = ("\nUsage: drop_lines.pl --file1=<string> --file2=<string>
required:
\t\tfile1	Superset file
\t\tfile2	Set of lines to be dropped in the first file

All lines in file1 not being present in file2 will be parsed out and 
therefore create a distinct set to file2. Lines are defined by their 
first column. Output is written to STDOUT.\n
");
 
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "file1=s","file2=s");

	die("Please specify files\n") unless $CMD{file1};
	die("Please specify files\n") unless $CMD{file2};
  
	open FILE1, $CMD{file1} or die "Cannot open ".$CMD{file1}." for reading.\n";	
	open FILE2, $CMD{file2} or die "Cannot open ".$CMD{file2}." for reading.\n";
	
}
	      

