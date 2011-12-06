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
#  Module: Parser::FASTQ::fq_2_fa.pl
#  Purpose:
#  In:
#  Out:
#


# Convert FASTQ to FASTA
# written by Korbinian Schneeberger, Stephan Ossowski

use Getopt::Long;

my %CMD;
my $fastq = "";

GetCom();


open FASTQ, $fastq or die "Cannot open fastq file\n";

my $c = 1;

while (my $line = <FASTQ>) {
	if(substr($line, 0, 1) eq "@") {
		my @a = split " ", $line;
		#print ">".$a[0]."\n";
		print ">".$c."\n";
		my $seq = <FASTQ>;
		chomp($seq);
		print "$seq\n";
		$seq = <FASTQ>;
		$seq = <FASTQ>;
		$c++;
	}
	else { print "file format not correct!\n"; exit(0); }
}




exit(0);

sub GetCom{

  my @usage = ("Usage: $0\n

required:
--fastq\tfastq formatted file

Will be converted into a fastq file.
	");

	die @usage if ($ARGV[0] eq "");
	GetOptions(\%CMD, "fastq=s");

	die("Please specify fastq file\n") unless defined($CMD{fastq});
  
	$fastq = $CMD{fastq};


}
	      

