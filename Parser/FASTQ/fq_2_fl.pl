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
#  Module: Parser::FASTQ::fq_2_fl.pl
#  Purpose:
#  In:
#  Out:
#


# Convert FASTQ to Shore sequence read format
# Written by Korbinian Schneeberger, Stephan Ossowski


use Getopt::Long;

my %CMD;
my $fastq = "";
my $length = "";

GetCom();


### prepare fake qcal and chas
my $qcal = ""; 
my $chas = "";
for(my $i = 0; $i < $length; $i++) {
	$qcal .= "Z";
	$chas .= "Z";
}

open FASTQ, $fastq or die "Cannot open fastq file\n";

while (<FASTQ>) {
	chomp;
	my $line = $_;

	if($line =~ /\@/) {
		my @members = split(/_/, $line);

		for (my $i = 3; $i < 6; $i++) {
			$members[$i] =~ s/-//g;
			if( length($members[$i]) == 1 ) { $members[$i] = "000" . $members[$i]; }
			elsif( length($members[$i]) == 2 ) { $members[$i] = "00" . $members[$i]; }
			elsif( length($members[$i]) == 3 ) { $members[$i] = "0" . $members[$i]; }
		}
		my $id = $members[2] . $members[3] . $members[4] . $members[5];

		my $seq = <FASTQ>;
		chomp($seq);
		print "$id\t$seq\t";
	}
	elsif($line =~ /^+/) {
		my $qval = <FASTQ>;
		chomp($qval);
		print "$qval\t$qcal\t$chas\n";
	}
	else { print "file format not correct!\n"; exit(0); }
}




exit(0);

sub GetCom{

  my @usage = ("Usage: $0\n

required:
--fastq\tfastq formatted file
--length\tlength of short reads


Will be converted into a fastq file.
	");

	die @usage if ($ARGV[0] eq "");
	GetOptions(\%CMD, "fastq=s", "length=s");

	die("Please specify fastq file\n") unless defined($CMD{fastq});
	die("Please specify length\n") unless defined($CMD{length});
  
	$fastq = $CMD{fastq};
	$length = $CMD{length};


}
	      

