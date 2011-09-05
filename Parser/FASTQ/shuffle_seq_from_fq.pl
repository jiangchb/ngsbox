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
#  Module: Parser::FASTQ::shuffle_seq_from_fq.pl
#  Purpose:
#  In:
#  Out:
#

#written by korbinian schneeberger and stephan ossowski

use Getopt::Long;

my %CMD;
my $num;
my $seq;
my $pe;

GetCom();

my $entries = 4;
if( $pe == 1) {
	$entries = 8;
}

my @fastq_entry;
while(<FILE>) {
	$fastq_entry[0] = $_;

	for(my $i = 1; $i < $entries; $i++) {
		my $line = <FILE>;
		$fastq_entry[$i] = $line;
	}

	my $rand = rand();
	if ($rand <= ($num / $seq)) {
      		$num--;
		for(my $i = 0; $i < $entries; $i++) {
      			print $fastq_entry[$i];
		}
	}
	$seq--;
}


sub GetCom{

  my @usage = ("Usage: $0 --file=<fq file> --num=<number of entries> 

required:
--file\t\tFq file to be parsed
--num\t\tNumber of sequences randomly taken out of the fq file
--perc\t\tPercent of sequences randomly taken out of the fq file (e.g. --perc 50)
--pe\t\tPaired end
\n"); 
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "file=s","num=s","perc=s","pe");

	die("Please specify a fl file\n") unless $CMD{file};
	die("Please specify number of entries\n") unless ($CMD{num} || $CMD{perc});

	open FILE, $CMD{file} or die "Cannot open file\n";	
	
	my $out = `wc -l $CMD{file}`;
	my @a = split " ", $out;
	$seq = $a[0];

	if($CMD{num}) {
		$num = $CMD{num};
	}
	elsif($CMD{perc}) {
		$num = int($seq * $CMD{perc} / 100);
	}

	if($CMD{pe}) {
		$pe = 1;
	}
	else {
		$pe = 0;
	}
}
