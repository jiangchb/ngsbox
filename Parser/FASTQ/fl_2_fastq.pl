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
#  Module: Parser::FASTQ::fl_2_fastq.pl
#  Purpose:
#  In:
#  Out:
#


### Convert Shore read files into fastq format
### written by Korbinian Schneeberger, Stephan Ossowski

use Getopt::Long;

my %CMD;
my $fl = "";

GetCom();

open FL, $fl or die "Cannot open fl file\n";

while (my $line = <FL>) {
	my @a = split " ", $line;

	print "@".$a[0]." ".$a[2]."\n";
	print $a[1]."\n";
	print "+\n";
	print $a[3]. "\n";
}



exit(0);

sub GetCom{

  my @usage = ("Usage: $0\n

required:
--fl\tfl formatted file

Will be converted into a fastq file.
	");

	die @usage if ($ARGV[0] eq "");
	GetOptions(\%CMD, "fl=s");

	die("Please specify fl file\n") unless defined($CMD{fl});
  
	$fl = $CMD{fl};
}
	      

