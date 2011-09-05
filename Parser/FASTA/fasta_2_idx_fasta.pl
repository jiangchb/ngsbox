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
#  Module: Parser::FASTA::fasta_2_idx_fasta.pl
#  Purpose:
#  In:
#  Out:
#

use Getopt::Long;

my %CMD;
my $file;
my $out;
my $idx;

GetCom();

open FILE, $file or die "Need a fasta file\n";
open OUT, ">".$out or die "Need a fasta file\n";
open IDX, ">".$idx or die "Need a fasta file\n";

my $count = 1;
while (<FILE>) {
	chomp();
	if (substr($_, 0, 1) eq ">") {
		my @a = split " ", $_;
		print OUT ">$count\n";
		#print IDX $count, "\t", substr($a[0], 1, length($a[0])-1), "\n";
		print IDX $count, "\t", substr($_, 1, length($_)-1), "\n";
		$count++;
	} else {
		print OUT $_, "\n";
	}
}


sub GetCom{

        my @usage = ("\nUsage: $0 --file=file --out=filename --idx=filename
                --file\treference genome
		--out\treference genome with index in fasta headers (used for mapping)
		--idx\t1 to 1 look-up for reference file (original fasta entry to idx fasta entry)

                description:

                \n\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "file=s", "out=s", "idx=s");

        die("Please specify file\n")    if not defined $CMD{file};
	die("Please specify out\n")    if not defined $CMD{out};
	die("Please specify idx\n")    if not defined $CMD{idx};

        $file = $CMD{file};
	$out = $CMD{out};
	$idx = $CMD{idx};

}



