#!/usr/bin/perl

# --------------------------------------------------------------------------
# Sliding window analysis of call-ability measured by percent reference call
# Written by Stephan Ossowski
# --------------------------------------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use FindBin;

### Command line parameters -------------------------------------------------------------------
my $ref      = "";
my $chrsizes = "";
my $win      = 200000;
my $win_step = 10000;

my %CMD;
GetCom();

my $slw_file = "$ref.slw";
my $pdf_file = "$ref.pdf";


### Sliding window ref-call analysis ------------------------------------------------------
my $chr = "NA";
my $genome_pos = 1;
my $front_start = $win_step - $win;
my %ref_sums = ();

open SLW, ">$slw_file" or die "Cannot open output file.\n";
open REF, $ref or die "Cannot open ref file.\n";
while( <REF> ) {
	chomp;
	my @a = split("\t", $_);

	### New chromosome reached
	if( $chr ne $a[1] ) {

		# First chromosome: Initialize
		if($chr eq "NA") {
			$chr = $a[1];
			$genome_pos = 1;
			$front_start = -$win;
		}

		# Later chromosomes
		else {
			# Print last sliding window of chromosome
			if($front_start >= 0) {
				if(! exists $ref_sums{$front_start}) { $ref_sums{$front_start} = 0; }
				$ref_sums{$front_start} /= $win;
				my $start = $front_start + 1;
				my $end   = $front_start + $win;
				print SLW "$chr\t$start\t$end\t" . $ref_sums{$front_start} . "\n";
			}

			# Reset
			$chr = $a[1];
			$genome_pos = 1;
			$front_start = $win_step - $win;
			%ref_sums = ();

			# Small Alyr fix: only use first 8 chr"
			my ($junk, $scaffold) = split("_", $chr);
			if($scaffold > 8) { last; }
		}
	}

	### Walk along the chromosome
	while( $genome_pos < $a[2] ) {
		$genome_pos++;

		# Sliding window end reached: Store and remove window
		if(  ($genome_pos % $win_step) == 0 ) {
			if($front_start >= 0) {
				if(! exists $ref_sums{$front_start}) { $ref_sums{$front_start} = 0; }
				$ref_sums{$front_start} /= $win;
				my $start = $front_start + 1;
				my $end   = $front_start + $win;
				print SLW "$chr\t$start\t$end\t" . $ref_sums{$front_start} . "\n";
			}
			
			delete $ref_sums{$front_start};
			$front_start += $win_step;
		}
	}

	### Update ref_sums
	for( my $i = $front_start; $i < $front_start + $win; $i += $win_step) {
		$ref_sums{$i}++;
	}
}
close(REF);



### Call R for plotting -----------------------------------------------------------------------
my $cmd = "R --slave --vanilla --args $chrsizes $pdf_file $slw_file $win_step $win < $FindBin::Bin/callability.R";
print STDERR $cmd, "\n";
system($cmd);

exit(0);



### Read command line parameters --------------------------------------------------------------
sub GetCom {
  my @usage = ("\nUsage: $0

--refseq     STRING     Reference calls in Shore format (currently reference.txt)
--chrsizes   STRING     Chromosome sizes file
--winsize    INT        Sliding window size
--winstep    INT        Sliding window step size
\n");


	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "refseq=s", "chrsizes=s", "winsize=s", "winstep=s");

	die("Please specify refseq file\n") unless defined($CMD{refseq});
	die("Please specify chr sizes file\n") unless defined($CMD{chrsizes});
	die("Please specify window size\n") unless defined($CMD{winsize});
	die("Please specify step size\n") unless defined($CMD{winstep});

	$ref  = $CMD{refseq};
	$chrsizes = $CMD{chrsizes};
	$win = $CMD{winsize};
	$win_step = $CMD{winstep};
}

