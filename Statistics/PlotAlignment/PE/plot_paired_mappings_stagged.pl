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
#  Module: Statistics::PlotAlignment::PE::plot_paired_mappings_stagged.pl
#  Purpose:
#  In:
#  Out:
#

####################################################################################
#Author 	Korbinian Schneeberger 
#Date 		05/15/07
#Version	0.1
#Input		map.list
#Function	Plots a user specified mapping regions (Yieha)
####################################################################################

use Getopt::Long;
use Cwd;

### Get command line options
my %CMD = ();

my $file;
my $idx;
my $chr;
my $begin;
my $end;

my $suppress_singletons = 0;
my $suppress_repeats    = 0;
my $nobackbone          = 0;

my $maxinsertsize = 0;
my $seqlength     = 35;
my $gap_size      = 1;
my $database      = "solexa";
my $format        = "pdf";

my $dir = getcwd;
GetCom();


### Initialize other variables
my $curr_pos;

my $min_x = -1;
my $max_x = -1;
my $max_y =  1;

my $min_pos = $begin - $maxinsertsize;
my $max_pos = $end + $maxinsertsize;
$begin -= $maxinsertsize;	#DEBUG
$end += $maxinsertsize;		#DEBUG
$min_pos = 1 if $min_pos < 1;

my %colors = (
	"u0" => "grey", "r0" => "grey", "u1" => "grey", "r1" => "grey", "u2" => "grey", "r2" => "grey",
	"u3" => "limegreen", "r3" => "darkgreen" , "u4" => "orange", "r4" => "darkorange3", "u5" => "red", "r5" => "darkred",
	"u6" => "limegreen", "r6" => "darkgreen" , "u7" => "orange", "r7" => "darkorange3", "u8" => "red", "r8" => "darkred",
	#"u9" => "limegreen", "r9" => "darkgreen" , "u10" => "orange", "r10" => "darkorange3", "u11" => "red", "r11" => "darkred",
	#"u12" => "limegreen", "r12" => "darkgreen" , "u13" => "orange", "r13" => "darkorange3", "u14" => "red", "r14" => "darkred"
	"u9" => "lightblue", "r9" => "darkblue" , "u10" => "gold", "r10" => "gold4", "u11" => "deeppink", "r11" => "deeppink3",
	"u12" => "lightblue", "r12" => "darkblue" , "u13" => "gold", "r13" => "gold4", "u14" => "deeppink", "r14" => "deeppink3"
);

my %paired_reads = ();


### Open input and output files
open FILE, $file or die "Cannot open $file\n";
open IDX, $idx or die "Cannot open $idx\n";
open MAPPING, "> $dir/tmpmap_paired_4_plotting.txt" or die "Cannot open tmp file\n";
open BRIDGE, "> $dir/tmpmap_bridge_4_plotting.txt" or die "Cannot open tmp file\n";
open MATRIX, "> $dir/tmpmap_matrix_4_plotting.txt" or die "Cannot open tmp file\n";


### Read in index file and jump FILE to right position
my $jump = 0;
while (my $l = <IDX> ) {
	my @a = split " ", $l;
	if ($a[0] > $chr or ($a[0] == $chr and $a[1] >= $min_pos)) {
		last;
	} else {
		$jump = $a[2];
	}
}
seek(FILE, $jump, 0);

print STDERR "JUMPED TO: $jump\n";


### Find first mapped read after (or at) starting posion
while (<FILE>) {
	my @a = split " ", $_;

	### Check if starting pos is reached
	if ($a[0] == $chr and $a[1] >= $min_pos and ($suppress_repeats == 0 or $a[6] == 1)) {

		### Add new read ID to read-hash if not exists
		if(! exists $paired_reads{$a[3]}) {
			my %tmp = ();
			$paired_reads{$a[3]} = \%tmp;
		}

		### Set min and max position of read or read pair
		if( (! exists $paired_reads{$a[3]}{min}) || ($a[1] < $paired_reads{$a[3]}{min}) ) {
			$paired_reads{$a[3]}{min} = $a[1];
		}
		if( (! exists $paired_reads{$a[3]}{max}) || ( ($a[1] + $a[7]) > $paired_reads{$a[3]}{max}) ) {
			$paired_reads{$a[3]}{max} = $a[1] + $a[7];
		}

		### Store read data
		if(	$a[9] == 0 || $a[9] == 1  ||			# single end or unresolved 
			$a[9] == 3 || $a[9] == 4  || $a[9] == 5 ||	# paired end
			$a[9] == 9 || $a[9] == 10 || $a[9] == 11	# mate pairs
		) {
			$paired_reads{$a[3]}{1} = \@a;
		}
		elsif(	$a[9] == 2  ||					# unresolved
			$a[9] == 6  || $a[9] == 7  || $a[9] == 8 ||	# paired end
			$a[9] == 12 || $a[9] == 13 || $a[9] == 14	# mate pairs	
		) {
			$paired_reads{$a[3]}{2} = \@a;
		}

		$curr_pos = $a[1];
		$min_x = $a[1];
		last;
	}
}

print STDERR "FOUND FIRST ENTRY: $curr_pos\n";


### Get all reads from region
my $line_counter = 0;
while ($curr_pos <= $max_pos) {
        my $line = <FILE>;

$line_counter++;

        if (not defined($line)) {
		last;
	}
        else {
                my @a = split " ", $line;

		### Check if end of region is reached
		if ($a[0] == $chr and $a[1] <= $max_pos) {

			if ($suppress_repeats == 0 or $a[6] == 1) {

				### Add new read ID to read-hash if not exists
				if(! exists $paired_reads{$a[3]}) {
					my %tmp = ();
					$paired_reads{$a[3]} = \%tmp;
				}

				### Set min and max position of read or read pair
				if( (! exists $paired_reads{$a[3]}{min}) || ($a[1] < $paired_reads{$a[3]}{min}) ) {
					$paired_reads{$a[3]}{min} = $a[1];
				}
				if( (! exists $paired_reads{$a[3]}{max}) || ( ($a[1] + $a[7]) > $paired_reads{$a[3]}{max}) ) {
					$paired_reads{$a[3]}{max} = $a[1] + $a[7];
				}

				### Store read data
				if(	$a[9] == 0 || $a[9] == 1  ||			# single end or unresolved
					$a[9] == 3 || $a[9] == 4  || $a[9] == 5 ||	# paired end
					$a[9] == 9 || $a[9] == 10 || $a[9] == 11	# mate pairs
				) {
					$paired_reads{$a[3]}{1} = \@a;
				}
				elsif(	$a[9] == 2  ||					# unresolved
					$a[9] == 6  || $a[9] == 7  || $a[9] == 8 ||	# paired end
					$a[9] == 12 || $a[9] == 13 || $a[9] == 14	# mate pairs
				) {
					$paired_reads{$a[3]}{2} = \@a;
				}
			}
		}

                $curr_pos = $a[1];
		$max_x = $a[1] + $a[7];
        }
}

print STDERR "Number of reads:". $line_counter."\n";
print STDERR "FINISHED READING ENTRIES: $curr_pos\n";


### Create plotting matrix
my %read_matrix   = ();
my %plot_matrix   = ();
my %bridge_matrix = ();
my %bridge_color  = ();
my %y_blocked     = ();
$y_blocked{1} = 0;

foreach my $read_id ( sort {$paired_reads{$a}{min} <=> $paired_reads{$b}{min}} keys %paired_reads ) {

	### Exclude singleton reads if option is selected.
	if(     ( exists $paired_reads{$read_id}{1} && exists $paired_reads{$read_id}{2} ) ||
		($suppress_singletons == 0) 
	) {
	
		### Check if there is a free column for the new read in the matrix
		my $use_y = -1;
		for( my $y = 1; $y <= $max_y; $y += 5 ) {
			if( $y_blocked{$y} < $paired_reads{$read_id}{min} ) {
				$use_y = $y;
				$y_blocked{$use_y} = $paired_reads{$read_id}{max} + 15;
				last;
			}
		}

		### No free column found in matrix, extend matrix
		if($use_y == -1) {
			$max_y += 5;
			$use_y = $max_y;
			$y_blocked{$use_y} = $paired_reads{$read_id}{max} + 15;
		}

		### Store reads and bridge between paired reads. Exclude singleton reads if option is selected.
		my $rep = "";

		# Read 1
		if( exists $paired_reads{$read_id}{1} ) {
			### Check if read is repetitive
			$rep = "u";
			if($paired_reads{$read_id}{1}[6] > 1) {
				$rep = "r";
			}

			### Print read wise matrix
			$read_matrix{$paired_reads{$read_id}{1}[1] . "," . $use_y} = $paired_reads{$read_id}{1}[4] . $rep . $paired_reads{$read_id}{1}[9];

			### Print position wise matrix : Currently commented for speed and testing reasons
			#for( my $x = $paired_reads{$read_id}{1}[1]; $x < $paired_reads{$read_id}{1}[1] + $paired_reads{$read_id}{1}[7]; $x++ ) {
			#	$plot_matrix{"$x,$use_y"} = $paired_reads{$read_id}{1}[4] . $rep . $paired_reads{$read_id}{1}[9];
			#}
		}
		# Read 2
		if( exists $paired_reads{$read_id}{2} ) {
			### Check if read is repetitive
			$rep = "u";
			if($paired_reads{$read_id}{2}[6] > 1) {
				$rep = "r";
			}

			### Print read wise matrix
			$read_matrix{$paired_reads{$read_id}{2}[1] . "," . $use_y} = $paired_reads{$read_id}{2}[4] . $rep . $paired_reads{$read_id}{2}[9];

			### Print position wise matrix : Currently commented for speed and testing reasons
			#for( my $x = $paired_reads{$read_id}{2}[1]; $x < $paired_reads{$read_id}{2}[1] + $paired_reads{$read_id}{2}[7]; $x++ ) {
			#	$plot_matrix{"$x,$use_y"} = $paired_reads{$read_id}{2}[4] . $rep . $paired_reads{$read_id}{2}[9];
			#}
		}

		# Bridge between paired reads
		if( exists $paired_reads{$read_id}{1} && exists $paired_reads{$read_id}{2} ) {
			my $start = 0;
			my $end = 0;
			if($paired_reads{$read_id}{1}[4] eq "D") {
				$start = $paired_reads{$read_id}{1}[1] + $seqlength; #$paired_reads{$read_id}{1}[7]; <- use this for real matrix, seqlength is only an approx
				$end   = $paired_reads{$read_id}{2}[1];
			}
			else {
				$start = $paired_reads{$read_id}{2}[1] + $seqlength; #$paired_reads{$read_id}{2}[7]; <- use this for real matrix, seqlength is only an approx
				$end   = $paired_reads{$read_id}{1}[1];
			}
	
			$bridge_matrix{"$start,$use_y"} = "$end,$use_y";
			$bridge_color{"$start,$use_y"} = $paired_reads{$read_id}{1}[4] . $rep . $paired_reads{$read_id}{1}[9];
		}
	}
}

print STDERR "FINSISHED PLOTTING MATRIX\n";


### Print plot_matrix : Currently commented for speed and testing reasons
#for (my $y = 1; $y <= $max_y; $y += 5) {
#	for (my $x = $min_x; $x <= $max_x; $x++) {
#		if(exists $plot_matrix{"$x,$y"}) {
#			print MATRIX "1";
#		}
#		else {
#			print MATRIX "0";
#		}
#	}
#	print MATRIX "\n";
#}

### Print read_matrix
foreach my $key ( keys %read_matrix) {
	my ($x, $y) = split(",", $key);
	print MAPPING "$x\t$y\t" . substr($read_matrix{$key}, 0, 1) . "\t" . $colors{ substr($read_matrix{$key}, 1) } . "\n";
}

### Print bridge matrix
foreach my $key ( keys %bridge_matrix) {
	my ($x1, $y1) = split(",", $key);
	my ($x2, $y2) = split(",", $bridge_matrix{$key});
	print BRIDGE "$x1\t$y1\t$x2\t$y2\t" . $colors{ substr($bridge_color{$key}, 1) } . "\n";
}

close MAPPING;
close BRIDGE;
close MATRIX;

print STDERR "FINISHED\n";

### Plot using R
my $r_cmd = "R --slave --vanilla --args '$dir/tmpmap_paired_4_plotting.txt' '$dir/tmpmap_bridge_4_plotting.txt' $format $seqlength $chr $begin $end < ".$ENV{PGSP}."/Plot/PE/plot_mapping.R";
print STDOUT $r_cmd, "\n";
system("$r_cmd");
system("rm $dir/tmpmap_paired_4_plotting.txt $dir/tmpmap_bridge_4_plotting.txt $dir/tmpmap_matrix_4_plotting.txt");

exit(0);


sub GetCom {

  my @usage = ("\nUsage: $0 

required:
--file\t\tmap.list file produced by SHORE
--chr\t\tChromosome
--begin\t\tPlot region begin position
--end\t\tPlot region end position

optional:
-suppress_singletons\t\tDo not plot singleton reads or reads missing their partner within the plotting interval
-suppress_repeats\t\tDo not plot repetitive reads
-nobackbone\tNo backbone available (currently always on until database and seq_ref table are implemented)

--maxinsertsize\t\tPaired end insert size
--idx\t\tindex of map.list file (default: map.list filename concatenated with .idx)
--database\t\tMySQL database containing reference sequence table seq_ref, not yet implemented (default: solexa)
--format\t\tOutput format for plot file (png | pdf) (default: pdf)
\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "file=s", "chr=s", "begin=s", "end=s", "maxinsertsize=s", "idx=s", "database=s", "format=s", "nobackbone", "suppress_singletons", "suppress_repeats");

	die("Please specify an input file\n") unless $CMD{file};
	die("Please specify chr\n") unless $CMD{chr};
	die("Please specify begin\n") unless $CMD{begin};
	die("Please specify end\n") unless $CMD{end};

	$file = $CMD{file};
	$chr = $CMD{chr};
	$begin = $CMD{begin};
        $end = $CMD{end};

	if (defined($CMD{maxinsertsize})) {
		$maxinsertsize = $CMD{maxinsertsize};
	}
	if (defined($CMD{idx})) {
		$idx = $CMD{idx};
	} else {
		$idx = $file.".idx";
	}
	if (defined($CMD{database})) {
		$database = $CMD{database};
	}
	if (defined($CMD{format})) {
		$format = $CMD{format};
	}
	if (defined($CMD{nobackbone})) {
                $nobackbone = 1;
        }
	if (defined($CMD{suppress_singletons})) {
		$suppress_singletons = 1;
	}
	if (defined($CMD{suppress_repeats})) {
		$suppress_repeats = 1;
	}

	if ($format ne "png" and $format ne "pdf") {
		die "file format must either pdf or png\n";
	}

	if ($end < $begin) {
		die("end is smaller than begin... $end < $begin\n");
	}

	return(0);
}
