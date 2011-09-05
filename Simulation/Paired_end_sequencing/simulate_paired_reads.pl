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
#  Module: Simulation::Paired_end_sequencing::simulate_paired_reads.pl
#  Purpose:
#  In:
#  Out:
#


use Getopt::Long;

my %CMD;
my $file;
my $num;
my $dist;
my $readlen;
my @seq;

GetCom();

# Read in insert dist
open DIST, $dist or die "Cannot open dist file:",$dist,"\n";
my %ID = ();
my $dist_min = 10000000;
my $dist_max = 0;
my $num_dist = 0;
while (<DIST>) {
	my @a = split " ";
	$ID{$a[0]} = $a[1];
	$dist_min = $a[0] if $a[0] < $dist_min;
	$dist_max = $a[0] if $a[0] > $dist_max;
	$num_dist += $a[1];
}

# Read in fasta file
open FILE, $file or die "Cannot open fasta file\n";
my $seq_num = 0;
my $curr_seq = "";
while(<FILE>) {
	chomp($_);
	if (substr($_, 0, 1)  eq ">") {
		if ($curr_seq ne "") {
			push @seq, $curr_seq;
			$seq_num++;
		}		
		$curr_seq= "";
	} else {
		$curr_seq .= $_;
	}
}
push @seq, $curr_seq;


# Simulate paired read sequencing
open READ1, ">read1_ins.fl";
open READ2, ">read2_ins.fl";
open FASTA, ">reads.fl";
my $id = 9000200000000;
for (my $i = 0; $i<$num; $i++) {
	my $s = rand($seq_num); # select chromosome
	my $length = length($seq[$s]); # get length of chromosome
	if ($length >= $dist_max) {
		# keep distance to the end of the read
		my $pos1 = rand($length - $dist_max);
		my $read1 = substr($seq[$s], $pos1, $readlen);	
		
		# get pos for paired rea
		my $pos2 = get_paired_pos($pos1);
		my $read2 = substr($seq[$s], $pos2-$readlen, $readlen);
		my $read2_revcomp = rev_comp($read2);

		#print $read1, "\t", $read2_revcomp, "\n";
		$id++;
		my $flag = int(rand(2));
		if ($flag == 1) {
			print READ1 $id, "\t", $read1, "\t", 1, "\tQUAL1\tQUAL2\tQUAL3\n";
                        print READ2 $id, "\t", $read2_revcomp, "\t", 2, "\tQUAL1\tQUAL2\tQUAL3\n";
			print FASTA ">", $id, "\n", $read1, "\n>$id\n",  $read2_revcomp, "\n";
		}
		else {
			print READ1 $id, "\t", $read2_revcomp, "\t", 1, "\tQUAL1\tQUAL2\tQUAL3\n";
			print READ2 $id, "\t", $read1, "\t", 2, "\tQUAL1\tQUAL2\tQUAL3\n";
			print FASTA ">", $id, "\n", $read2_revcomp, "\n>$id\n",  $read1, "\n";
		}
	} else {
		print STDERR ("Found seq too short for paired read sequencing -- refusing output");
	}

}

sub get_paired_pos {
	my ($pos1) = @_;

	my $rand = rand($num_dist)+1;

	my $pos2 = 0;
	my $i = $dist_min - 1;

	while ($rand > 0) {
		$i++;
		$pos2 = $i;
		if (defined($ID{$i})) {
			$rand -= $ID{$i};
		}
		if ($i > $dist_max) {
			last;
		}
	}

	return $pos2 + $pos1;
}


sub rev_comp {
        my ($seq) = @_;
        my $new = reverse($seq);
        $new =~ tr/acgtACGT/tgcaTGCA/;

        return $new;
}



sub GetCom{
  my @usage = ("Usage: $0 --fasta=<file> --num=<int> --dist=<int> --dev=<int> --readlen=<int>

required:
--fasta\t\tDefine fasta file
--num\t\tDefine number of read pairs
--dist\t\tinsert_dist.txt file
--readlen\tDefined sequencing length
optional:
");
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "fasta=s", "num=s", "dist=s",  "readlen=s");

	die("Please specify a fasta file\n") unless $CMD{fasta};
	die("Please specify number of read pairs\n") unless $CMD{num};
	die("Please specify insert size of fragment\n") unless $CMD{dist};
	die("Please specify sequencing length\n") unless $CMD{readlen};
  
	$file = $CMD{fasta};
	$num = $CMD{num};
	$dist = $CMD{dist};
	$readlen = $CMD{readlen};

}
	      

