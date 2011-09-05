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
#  Module: Statistics::PlotAlignment::PE::plot_paired_mappings.pl
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

my %CMD = ();
my $dir = getcwd;
my $idx;
my $file;
my $begin;
my $end;
my $chr = 1;
my $database = "solexa";
my $format = "pdf";
my $maxinsertsize;
my $seqlength = 35;
my $gap_size = 1;
my $nobackbone;

GetCom();

open FILE, $file or die "Cannot open $file\n";
open IDX, $idx or die "Cannot open $idx\n";

open MAPPING, "> $dir/tmpmap_paired_4_plotting.txt" or die "Cannot open tmp file\n";
open BRIDGE, "> $dir/tmpmap_bridge_4_plotting.txt" or die "Cannot open tmp file\n";

my $min_pos = $begin - $maxinsertsize;
my $max_pos = $end + $maxinsertsize;
$min_pos = 1 if $min_pos < 1;
my $curr_pos;

my $paired_single = "limegreen";
my $paired_rep = "darkgreen";
my $unpaired_single = "red";
my $unpaired_rep = "darkred";
my $pairanomaly_single = "orange";
my $pairanomaly_rep = "darkorange3";

my @MAPPING_POS = ();

my @MAPPING_X = ();
my @MAPPING_Y = ();
my @MAPPING_DIR = ();
my @MAPPING_COL = ();

my %PAIR_X = (); # stores the ids at a given (key) position
my %PAIR_DIR = (); # stores strand by id 
my %PAIR_COL = (); # stores color by id

my %BRIDGE_X_1 = ();
my %BRIDGE_Y_1 = ();
my %BRIDGE_X_2 = ();
my %BRIDGE_Y_2 = ();
my %BRIDGE_COL = ();

############# Read in index file and jump FILE to right position ############
my $jump = 0;
SETTING: while (my $l = <IDX> ) {
	my @a = split " ", $l;
	if ($a[0] > $chr or ($a[0] == $chr and $a[1] >= $min_pos)) {
		last SETTING;
	} else {
		$jump = $a[2];
	}
}
seek(FILE, $jump, 0);

print STDERR "GOT IT $jump\n";

# set file to region of interest and read in first fragment
INIT: while (<FILE>) {
	my @a = split " ", $_;
	if ($a[0] == $chr and $a[1] >= $min_pos) {

		$curr_pos = $a[1];
		$MAPPING_POS[0] = ($curr_pos - $min_pos) + $seqlength + $gap_size;

		$MAPPING_X[0] = $a[1];
		$MAPPING_DIR[0] = $a[3];
		$MAPPING_Y[0] = 1;

		my $code = $a[7];
		my $pair_x = $a[9];
		my $pair_dir = $a[11];
		my $read_id = $a[5];
		my $read_dir = $a[3];

		if (substr($code, 0, 1) eq "P") {
 			if (substr($code, 1, 1) eq "U") {
				$MAPPING_COL[0] = $paired_single;
				set_bridge($read_id, $a[1], 1, $read_dir, $paired_single);
			}
			elsif (substr($code, 1, 1) eq "R") {
				$MAPPING_COL[0] = $paired_rep;
			}

			$PAIR_X{$pair_x} .= $read_id."#";
			$PAIR_DIR{$read_id} = $pair_dir;
			if (substr($code, 2, 1) eq "U") {
				$PAIR_COL{$read_id} = $paired_single;
			}
			elsif (substr($code, 2, 1) eq "R") {
				$PAIR_COL{$read_id} = $paired_rep;
			}

		}
		elsif (substr($code, 0, 1) eq "A") {
			if (substr($code, 1, 1) eq "U") {
                                $MAPPING_COL[0] = $pairanomaly_single;
				set_bridge($read_id, $a[1], 1, $read_dir, $pairanomaly_single);
                        }
                        elsif (substr($code, 1, 1) eq "R") {
                                $MAPPING_COL[0] = $pairanomaly_rep;
                        }
			

                        $PAIR_X{$pair_x} = $read_id."#";
                        $PAIR_DIR{$read_id} = $pair_dir;
                        if (substr($code, 2, 1) eq "U") {
                                $PAIR_COL{$read_id} = $pairanomaly_single;
                        }
                        elsif (substr($code, 2, 1) eq "R") {
                                $PAIR_COL{$read_id} = $pairanomaly_rep;
                        }

		}
		elsif (substr($code, 0, 1) eq "M") {
			if (substr($code, 1, 1) eq "U") {
                                $MAPPING_COL[0] = $unpaired_single;
                        }
                        elsif (substr($code, 1, 1) eq "R") {
                                $MAPPING_COL[0] = $unpaired_rep;
                        }
		}
			

		last INIT;
	}
}

print STDERR "FILE SET\n";

# Parse until the end of region of interest
my $distance = 0;
my $count = 0;
while ($curr_pos <= $max_pos) {
	$count++;

	my $line = <FILE>;
	my @a;
	my $new_curr_pos;

	if (not defined($line)) {
		$new_curr_pos = $max_pos + 1;
	}
	else {
		@a = split " ", $line;
		$new_curr_pos = $a[1];
	}

	### Empty paired mappings storage
	for (my $pos = $curr_pos; $pos <= $new_curr_pos; $pos++) {	
	        if (defined($PAIR_X{$pos})) {

			my @ids = split "#", $PAIR_X{$pos};

			$distance = $pos - $min_pos;
			## All pairs starting at this position
			foreach my $id (@ids) {
				push @MAPPING_X, $pos;
				push @MAPPING_COL, $PAIR_COL{$id};
				push @MAPPING_DIR, $PAIR_DIR{$id};

				SET: for (my $y = 0; $y<=@MAPPING_POS; $y++) {
					if (not defined($MAPPING_POS[$y]) or $MAPPING_POS[$y] <= $distance) {
						push @MAPPING_Y, $y + 1;
						$MAPPING_POS[$y] = $distance + $seqlength + $gap_size;
						
						if ($PAIR_COL{$id} eq $paired_single or $PAIR_COL{$id} eq $pairanomaly_single) {						
							set_bridge($id, $pos, $y+1, $PAIR_DIR{$id}, $PAIR_COL{$id});
						}
	
		        	                last SET;
        		        	}
	        		}
				delete $PAIR_COL{$id}; # delete $PAIR_DIR{$id};
			}
			delete $PAIR_X{$pos};
		}
        }

	$curr_pos = $new_curr_pos;
	$distance = $curr_pos - $min_pos;

	### Write out current line
	
	if ($curr_pos <= $max_pos) {

		if (not defined($PAIR_DIR{$a[5]})) { # Prevent the pair acts as first read again

			$curr_pos = $a[1];
			my $curr_y;

        	        push @MAPPING_X, $a[1];
                	push @MAPPING_DIR, $a[3];

			SET: for (my $y = 0; $y<=@MAPPING_POS; $y++) {		
				if (not defined($MAPPING_POS[$y]) or $MAPPING_POS[$y] <= $distance) {
					push @MAPPING_Y, $y + 1;
                        	        $MAPPING_POS[$y] = $distance + $seqlength + $gap_size;
					#$curr_y = $MAPPING_POS[$y];
					$curr_y = $y + 1;
                                	last SET;
				}
			}

			my $code = $a[7];
			my $pair_x = $a[9];
			my $pair_dir = $a[11];
	                my $read_id = $a[5];
        	        my $read_dir = $a[3];

			if (substr($code, 0, 1) eq "P") {
        	                if (substr($code, 1, 1) eq "U") {
                	                push @MAPPING_COL,  $paired_single;
					set_bridge($read_id, $curr_pos, $curr_y, $read_dir, $paired_single);
                        	}
	                        elsif (substr($code, 1, 1) eq "R") {
        	                        push @MAPPING_COL, $paired_rep;
                	        }


	                        $PAIR_X{$pair_x} = $read_id."#";
        	                $PAIR_DIR{$read_id} = $pair_dir;
                	        if (substr($code, 2, 1) eq "U") {
                        	        $PAIR_COL{$read_id} = $paired_single;
	                        }
        	                elsif (substr($code, 2, 1) eq "R") {
                	                $PAIR_COL{$read_id} = $paired_rep;
                        	}

	                }
        	        elsif (substr($code, 0, 1) eq "A") {
                	        if (substr($code, 1, 1) eq "U") {
                        	        push @MAPPING_COL, $pairanomaly_single;
					set_bridge($read_id, $curr_pos, $curr_y, $read_dir, $pairanomaly_single);
	                        }
        	                elsif (substr($code, 1, 1) eq "R") {
                	                push @MAPPING_COL, $pairanomaly_rep;
                        	}

	
        	                $PAIR_X{$pair_x} = $read_id."#";
                	        $PAIR_DIR{$read_id} =  $pair_dir;
                        	if (substr($code, 2, 1) eq "U") {
                                	$PAIR_COL{$read_id} = $pairanomaly_single;
	                        }
        	                elsif (substr($code, 2, 1) eq "R") {
                	                $PAIR_COL{$read_id} = $pairanomaly_rep;
                        	}

	                }
	                elsif (substr($code, 0, 1) eq "M") {
        	                if (substr($code, 1, 1) eq "U") {
                	                push @MAPPING_COL, $unpaired_single;
	                        }
        	                elsif (substr($code, 1, 1) eq "R") {
                	                push @MAPPING_COL, $unpaired_rep;
                        	}
	                }
		}
	}
}


#######################################################
#
# Print mapping and pairings
#
for (my $i = 0; $i<@MAPPING_X; $i++) {
	print MAPPING $MAPPING_X[$i], "\t", $MAPPING_Y[$i], "\t", $MAPPING_DIR[$i], "\t", $MAPPING_COL[$i],"\n";
}

my $count_keys = 0;
my $count_bridges = 0;
foreach my $id (keys %BRIDGE_X_1) {
	$count_keys++;
	if (defined($BRIDGE_X_1{$id}) and defined($BRIDGE_X_2{$id})) {
		$count_bridges++;
		print BRIDGE $BRIDGE_X_1{$id}, "\t", $BRIDGE_Y_1{$id}, "\t", $BRIDGE_X_2{$id}, "\t", $BRIDGE_Y_2{$id}, "\t", $BRIDGE_COL{$id}, "\n";
	}
}
print STDERR "Tried bridges: ", $count_keys, "\tBridges build: ", $count_bridges, "\n";

close MAPPING;
close BRIDGE;

# Plot using R
system("R --slave --vanilla --args '$dir/tmpmap_paired_4_plotting.txt' '$dir/tmpmap_bridge_4_plotting.txt' $format $seqlength $chr $begin $end < /ebio/abt6/korbinian/pgsp/Plot/PE/plot_mapping.R");
system("rm $dir/tmpmap_paired_4_plotting.txt $dir/tmpmap_bridge_4_plotting.txt");


sub set_bridge {
	my ($id, $x, $y, $dir, $col) = @_;

	if (not defined($BRIDGE_X_1{$id})) {
		if ($dir eq "D") {
			$BRIDGE_X_1{$id} = $x+$seqlength;
		}
		else {
			$BRIDGE_X_1{$id} = $x;	
		}
		$BRIDGE_Y_1{$id} = $y;
	}
	else {
		if ($dir eq "D") {
                        $BRIDGE_X_2{$id} = $x+$seqlength;
                }
                else {
                        $BRIDGE_X_2{$id} = $x;
                }
                $BRIDGE_Y_2{$id} = $y;
	}

	$BRIDGE_COL{$id} = $col;

	return(0);
}


sub GetCom {

  my @usage = ("\nUsage: $0 

required:
--file\t\tmap.list.pairs.sorted
--begin\t\tbegin position
--end\t\tend position
--chr\t\tMouse on Mars
--maxinsertsize\t\tPaired end insert

optional:
--idx\t\tindex of map.list.pairs.sorted (default: map.list filename concat with .idx)
--database\t\tdatabase (default: solexa)
--format\t\toutput file format (png | pdf) (default: pdf)
-nobackbone\tNo backbone available
\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "file=s", "chr=s", "begin=s", "end=s", "idx=s", "database=s", "format=s", "nobackbone", "maxinsertsize=s");

	die("Please specify an input file\n") unless $CMD{file};
	die("Please specify begin\n") unless $CMD{begin};
	die("Please specify end\n") unless $CMD{end};
	die("Please specify chr\n") unless $CMD{chr};
	die("Please specify insertsize\n") unless $CMD{maxinsertsize};

	$file = $CMD{file};
	$begin = $CMD{begin};
        $end = $CMD{end};
	$maxinsertsize = $CMD{maxinsertsize};

	if (defined($CMD{idx})) {
		$idx = $CMD{idx};
	} else {
		$idx = $file.".idx";
	}
	if (defined($CMD{chr})) {
		$chr = $CMD{chr};
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
	else {
		$nobackbone = 0;
	}
	

	if ($format ne "png" and $format ne "pdf") {
		die "file format must either pdf or png\n";
	}

	if ($end < $begin) {
		die("end is smaller than begin... $end < $begin\n");
	}

	return(0);
}



exit(0);
