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
#  Module: Analysis::SSR::ssr_accessibility_plot.pl
#  Purpose:
#  In:
#  Out:
#

# written by Korbinian Schneeberger and Stephan Ossowski


### Get user parameters
my $usage = "scan_ssr.pl <filename> <ssr length> <min. occ.> <rm periods>\nRepeats incl. \"N\" or \"n\" won't be reported!\n";

my $file       = shift 	or die $usage;		# reference file
my $refcall    = shift	or die $usage;		# quality reference call file
my $chrsize    = shift  or die $usage;		# chr size file
my $rep_size   = shift 	or die $usage;		# length of ssr: 1=homo, 2=di, 3=tri
my $min_rep    = shift 	or die $usage;		# minimum occurrencies (instances)
my $rm_periods = shift	or die $usage;		# flag: remove internal periods
my $type       = shift	or die $usage;		# Type = sequence of ssr or NA

my $cutoff = $rep_size * $min_rep or die $usage;	

open FILE, $file or die $usage;

### Check user parameters
if ($min_rep < 2) {
	print "Do you relly want to see 'repeats' with one instance?\n";
	exit(1);
}


### Variables and containers
my $seq = "";
my $id = "";

my %microsatelite = ();

my %microsat_occ = ();
my %microsat_avg_call = ();
my %microsat_hasgap =();


### Read in file
while (<FILE>) {
	chomp;
	if (substr($_, 0, 1) eq ">" ) {
		if ($seq ne "") {
			parse_seq($id, $rep_size);
			$seq = "";
		}
		my @a = split " ";
		$id = substr($a[0], 1, length($a[0])-1);
	} else {
		$seq .= $_;
	}
}

### Parse last chromosome
if ($seq ne "") {
	parse_seq($id, $rep_size);
}

print STDERR "Parsed " . scalar(keys %microsatelite) . " microsatellites\n\n";



### Parse reference callability
open CHRSIZE, $chrsize or die $usage;
while(<CHRSIZE>) {
	chomp;
	my ($chr, $size) = split("\t", $_);

	### Build refcall array
	my @chr_array = ();
	my $pos = 1;

	open RCALL, $refcall or die $usage;
	while( <RCALL> ) {
		chomp;
		my @e = split("\t", $_);
		if($chr == $e[0]) {
			
			while( $pos < $e[1] ) {
				$chr_array[$pos] = 0;
				$pos++;
			}
			$chr_array[$pos] = 1;
		}
	}
	while( $pos < $size) {
		$chr_array[$pos] = 0;
		$pos++;
	}

	### Check callability
	foreach my $instance ( sort {$a<=>$b} keys %microsatelite ) {
	
		# Get microsatellite data: id, start, end, base-count, instances, rep-type
		my @ms = split("-", $microsatelite{$instance});

		if( $ms[0] == $chr ) {
		
			# increase instance count
			$microsat_occ{$ms[4]}++;

			# check coverage along instance
			my $sum = 0;
			my $has_gap = 0;
			for( my $i = $ms[1]; $i <= $ms[2]; $i++) {
				$sum += $chr_array[$i];
				if( $chr_array[$i] ) {
					$has_gap = 1;
				}
			}

			# Check if instance has any callability gap
			if($has_gap == 1) { 
				$microsat_hasgap{$ms[4]}++;
			}
	
			# Calcualte average callability of instance
			my $avg_call = $sum / $ms[3];
			$microsat_avg_call{$ms[4]} += $avg_call;
		}
	
		if($ms[1] % 100 == 0) {
			print STDERR "$ms[0]\t$ms[1]\n";
		}
	}
}


### print results and run R
foreach my $instances (sort {$a<=>$b} keys %microsat_occ) {

	# Check if all hash entries are defined
	if(! exists $microsat_hasgap{$instances}) { $microsat_hasgap{$instances} = 0; }

	# Plot
	print $instances ."\t". $microsat_occ{$instances} ."\t". $microsat_hasgap{$instances} ."\t".
		$microsat_avg_call{$instances} / $microsat_occ{$instances} ."\n";
}

exit(0);



### Parse fasta entry for repeats
sub parse_seq {
	my ($id, $len) = @_;

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
						$rep = uc($rep);
						my $start = $i-$count+1;
						my $instances = int($count/$len);
						my $micro_string = $id . "-" . $start . "-" . $i . "-" . $count . "-" . $instances . "-" . $rep;

						if($type eq "NA" || $type eq $rep) {
							$microsatelite{"$id-$start"} = $micro_string;
							#print $id, "\t", $i-$count+1, "\t", $i, "\t", $count, "\t", int($count/$len), "\t", $rep, "\n";
						}
					}
				}
				else {
					$rep = uc($rep);
					my $start = $i-$count+1;
					my $instances = int($count/$len);
					my $micro_string = $id . "-" . $start . "-" . $i . "-" . $count . "-" . $instances . "-" . $rep;


					if($type eq "NA" || $type eq $rep) {
						$microsatelite{"$id-$start"} = $micro_string;
						#print $id, "\t", $i-$count+1, "\t", $i, "\t", $count, "\t", int($count/$len), "\t", $rep, "\n";
					}
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
