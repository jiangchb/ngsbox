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
#  Module: Analysis::SSR::scan_genome_4_ssr.pl
#  Purpose:
#  In:
#  Out:
#

# written by Korbinian Schneeberger


my $usage = "scan_ssr.pl <filename> <ssr length> <min. occ.> <rm periods>\nRepeats incl. \"N\" or \"n\" won't be reported!\n";
my $file       = shift 			or die $usage;
my $rep_size   = shift 			or die $usage;
my $min_rep    = shift 			or die $usage;
my $cutoff     = $rep_size * $min_rep 	or die $usage; # That is mononucleotids when scanning for tri-nucleotide repeats.
my $rm_periods = shift			or die $usage;

open FILE, $file 			or die $usage;

if ($min_rep < 2) {
	print "Do you relly want to see 'repeats' with one instance?\n";
	exit(1);
}


# read in file

my $seq = "";
my $id = "";
while (<FILE>) {
	chomp;
	if (substr($_, 0, 1) eq ">" ) {		
		if ($seq ne "") {
			parse_seq($id, $seq, $rep_size) ;
			$seq = "";
		}
		my @a = split " ";
		$id = substr($a[0], 1, length($a[0])-1);		
	} else {
		$seq .= $_;
	}
}
if ($seq ne "") {
	parse_seq($id, $seq, $rep_size) if $seq ne "";
	$seq = "";
}

# parse fasta entry for repeats

sub parse_seq {
	my ($id, $seq, $len) = @_;

	$seq =~ s/\s//g;
	$seq = uc($seq);

	my @window = ();
	
	init_window(\@window, $len);
	my $frame = 0;
	my $count = $len;


	for (my $i=0; $i<=length($seq)-$len; $i++) {
		if ($window[$frame] eq substr($seq, $i, 1)) {
			$count++;
		} else {
			my $rep = "";
			for (my $j=0; $j<@window; $j++) {$rep.=$window[$j];} # simulate join for speeding up
			if ($count >= $cutoff and $rep !~ "N") {
				if ($rm_periods == 1) {
					my $rm = 0;
					SEEK: for (my $k = 1; $k <= $rep_size/2; $k++) {
						my $periode_found = 1;	
						my $periode_inst = substr($rep, 0, $k);
						PER: for (my $p = 0; $p<length($rep); $p+=$k) {
							if (substr($rep, $p, $k) ne $periode_inst) {
								$periode_found = 0;
								last PER;	
							}
						}

						if ($periode_found == 1) {
							$rm = 1;
							last SEEK;
						}

					}
					if ($rm == 0) {
						print $id, "\t", $i-$count+1, "\t", $i, "\t", $count, "\t", int($count/$len), "\t", $rep, "\n";
					}
				}
				else {
					print $id, "\t", $i-$count+1, "\t", $i, "\t", $count, "\t", int($count/$len), "\t", $rep, "\n";
				}
			}
			$count = $len;
			$window[$frame] = substr($seq, $i, 1);
		}
		$frame++;
		$frame = $frame%$len;
	}	
}

sub init_window {
	my ($win_ref, $len) = @_;	

	for (my $i = 0; $i<$len; $i++) {
		@$win_ref[$i] = "";		
	}

}
