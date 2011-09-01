#! /usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use FindBin;
use lib $FindBin::Bin;

my $mapblt;

## Parse command lint
my %CMD;
GetCom();

my %DEL = ();
my %INS = ();
my %HDR = ();

open FILE, $mapblt or die "Cannot open file\n";
open INSINFO, ">map.blt.insertion.info" or die "cannot open file\n";
open HDRINFO, ">map.blt.hdr.info" or die "cannot open file\n";

while (my $line = <FILE>) {

	my @a = split " ", $line;

	my $target_chr = $a[13];
	my $ori = $a[8];

	my @block_sizes = split ",", $a[18];
	my @read_starts = split ",", $a[19];
	my @genome_starts = split ",", $a[20];

	# Only takes spanning alignments into account
	for (my $i = 0; $i < @block_sizes-1; $i++) {
		my $genome_start1 = 1 + $genome_starts[$i];
		my $genome_end1 = 1 + $genome_starts[$i] + $block_sizes[$i] - 1;
		my $genome_start2 = 1 + $genome_starts[$i+1];
                my $genome_end2 = 1 + $genome_starts[$i+1] + $block_sizes[$i+1] - 1;

		my $read_start1 = 1 + $read_starts[$i];
                my $read_end1 = 1 + $read_starts[$i] + $block_sizes[$i] - 1;
                my $read_start2 = 1 + $read_starts[$i+1];
                my $read_end2 = 1 + $read_starts[$i+1] + $block_sizes[$i+1] - 1;

		my $ins_length = $read_start2 - $read_end1 - 1;
		my $del_length = $genome_start2 - $genome_end1 - 1;

		if ($ins_length != 0 and $del_length != 0) {
			$HDR{$target_chr}{$genome_end1+1}{$genome_start2-1}{$ins_length}{$del_length}++;	
			print HDRINFO $target_chr, "\t", $genome_end1+1, "\t", $genome_start2-1, "\t", $del_length, "\t", $a[9], "\t", $read_end1+1, "\t", $read_start2-1, "\t", $ins_length, "\n";
		}
		elsif ($ins_length != 0) {
			$INS{$target_chr}{$genome_end1}{$ins_length}++;
			print INSINFO $target_chr, "\t", $genome_end1, "\t", $genome_end1+1, "\t0\t", $a[9], "\t", $read_end1+1, "\t", $read_start2-1, "\t", $ins_length, "\n";
		}
		else {
			$DEL{$target_chr}{$genome_end1+1}{$genome_start2-1}++;
		}

	}
}

close FILE;

open DELFILE, ">map.blt.deletion";
foreach my $chr (sort {$a <=> $b} keys %DEL) {
	foreach my $start (sort {$a <=> $b} keys %{$DEL{$chr}}) {
		foreach my $end (sort {$a <=> $b} keys %{$DEL{$chr}{$start}}) {
			print DELFILE $chr, "\t", $start, "\t", $end, "\t", $end-$start+1, "\t", "0", "\t", $DEL{$chr}{$start}{$end}, "\n";
		}
	}
}
close DELFILE;

open INSFILE, ">map.blt.insertion";
foreach my $chr (sort {$a <=> $b} keys %INS) {
        foreach my $start (sort {$a <=> $b} keys %{$INS{$chr}}) {
		foreach my $length (sort {$a <=> $b} keys %{$INS{$chr}{$start}}) {
	                print INSFILE $chr, "\t", $start, "\t", $start+1, "\t", "0", "\t", $length, "\t", $INS{$chr}{$start}{$length}, "\n";
		}
        }
}
close INSFILE;

open HDRFILE, ">map.blt.hdr";
foreach my $chr (sort {$a <=> $b} keys %HDR) {
        foreach my $start (sort {$a <=> $b} keys %{$HDR{$chr}}) {
        	foreach my $end (sort {$a <=> $b} keys %{$HDR{$chr}{$start}}) {
	                foreach my $ins_length (sort {$a <=> $b} keys %{$HDR{$chr}{$start}{$end}}) {
        	        	foreach my $del_length (sort {$a <=> $b} keys %{$HDR{$chr}{$start}{$end}{$ins_length}}) {
	        	                print HDRFILE $chr, "\t", $start, "\t", $end, "\t", $del_length, "\t", $ins_length, "\t", $HDR{$chr}{$start}{$end}{$ins_length}{$del_length}, "\n";
				}
			 }
                }
        }
}
close HDRFILE;





##########################################################################################
## parse anchor


sub GetCom {

        my @usage = ("$0 --mapblt file\n\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "mapblt=s");

        die("Please specify blat file\n") unless defined($CMD{mapblt});

        $mapblt = $CMD{mapblt};

}


