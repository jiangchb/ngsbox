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
#  Module: Parser::ML::shuffle_seq_from_map.pl
#  Purpose:
#  In:
#  Out:
#

#written by korbinian schneeberger

use Getopt::Long;

my %CMD;
my $num;
my $map;
my %ID = ();


GetCom();

open FILE, $map;
my $count_id = 0;

# Read in Read ids
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (!defined($ID{$a[3]})) {
		$count_id++;
	}
	$ID{$a[3]} = 0;
}
close FILE;

print STDERR "got read ids\n";

#  Print new map.list file
open FILE, $map;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (defined($ID{$a[3]}) {
		if ($ID{$a[3]} == 1) {
			print $line;
		}
	}
	else {
		my $rand = rand();
		if ($rand <= ($num / $count_id)) {
                	$num--;
        	        $ID{$a[3]} = 1;
			print $line;
	        }
		else {
			$ID{$a[3]} = 0;
		}
	        $count_id--;
	}
}


sub GetCom{

  my @usage = ("Usage: $0 --file=<map.list> --prob=<double> 

required:
--file\t\tmap.list file to be parsed
--prob\t\tProbability that a read is selected for subsampling
\n"); 
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "file=s","prob=s");

	die("Please specify a fl file\n") unless $CMD{file};
	die("Please specify number of entries\n") unless $CMD{num};
  
	$map = $CMD{file};
	$num = $CMD{num};

}
