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
#  Module: Statistics::PlotAlignment::SE::plot_mapping.pl
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
use DBI;

my $dbh;

my %CMD = ();
my $dir = getcwd;
my $idx;
my $file;
my $begin;
my $end;
my $chr = 1;
my $offset;
my $database = "solexa";
my $format = "pdf";
my $nobackbone;

GetCom();
if ($nobackbone == 0) {
	&connect_to_db();
}


open FILE, $file or die "Cannot open $file\n";
open IDX, $idx or die "Cannot open $idx\n";
open MAPPING, "> tmpmap_mapping_4_plotting.txt" or die "Cannot open tmp file\n";
open MM, "> tmpmap_mapping_mm.txt" or die "Cannot open tmp file\n";
open INSERTIONS, "> tmpmap_mapping_ins.txt" or die "Cannot open tmp file\n";
open SEQ, "> tmpmap_seq_ref.txt" or die "Cannot open seq ref file\n";
#open COVERAGE, "> mapping_cov_4_plotting.txt" or die "Cannot open mapping cov file for plotting\n";

my $start_pos = $begin - 120;
$start_pos = 1 if $start_pos < 1;
my $read_dist = 1;

my @MAPPING_POS = ();
my @MAPPING_X;
my @MAPPING_Y;
my @MAPPING_L;
my @MAPPING_COL;
my @MAPPING_COV;
my @MM_X;
my @MM_Y;
my @MM_C;
my @MM_PRB;
my @INS;
my @INS_PRB;
my $ins_max_x = -1;
my $ins_max_y = -1;

############# Read in index file an jump FILE to right position ############
my $jump = 0;
SETTING: while (my $l = <IDX> ) {
	my @a = split " ", $l;
	if ($a[0] > $chr or ($a[0] == $chr and $a[1] >= $start_pos)) {
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
	if ($a[0] == $chr and $a[1] >= $start_pos) {
		my $len = parse_seq($a[2], $a[1], 1, $a[10], $a[4]);
		$MAPPING_X[0] = $a[1];
		$MAPPING_Y[0] = 1;
		$MAPPING_POS[0] = $len + $read_dist;
		$MAPPING_L[0] = $len;
		$MAPPING_COL[0] = $a[6];
		last INIT;
	}
}

print STDERR "FILE SET\n";

# Parse until the end of region of interest
my $distance = 0;
my @set_array = ();

REG: while (<FILE>) {

	# Test for end of region
	my @a = split " ";
	if ($chr != $a[0] or $a[1] > $end) {
		last REG;
	}

	$distance = $a[1] - $start_pos;

	$MAPPING_COV[$distance]++;

	# Set x and y values for current read

############## OLD ########################
	SET: for (my $y = 0; $y<=@MAPPING_POS; $y++) {
		if (not defined($MAPPING_POS[$y]) or $MAPPING_POS[$y] <= $distance) {
			push @MAPPING_X, $a[1];
			push @MAPPING_Y, $y+1;
			my $seq_length  = parse_seq($a[2], $a[1], $y+1, $a[10], $a[4]);
			$MAPPING_POS[$y] = $distance + $seq_length + $read_dist;
			push @MAPPING_L, $seq_length;
			push @MAPPING_COL, $a[6];
			last SET;
		}
	}

############## NEW ########################
#	if ($distance >= ) {
#
#	}
#
#
#	Keep an array with the entries keep it order (intrinsic in the filling of the array) and get therefore the next
#	position, this can be done when holding the current coverage aswell, by counting the ins and outs of the array.
#	This could lead to the coverage plotting in addition.

}


#######################################################
#
# Print mapping, mismatches, insertions, coverage and backbone
#
for (my $i = 0; $i<@MAPPING_X; $i++) {
	if ($MAPPING_COL[$i] <= 3) {
		print MAPPING $MAPPING_X[$i], "\t", $MAPPING_Y[$i], "\t", $MAPPING_L[$i], "\t", $MAPPING_COL[$i],"\n";
	} else {
		print MAPPING $MAPPING_X[$i], "\t", $MAPPING_Y[$i], "\t", $MAPPING_L[$i], "\t", 4, "\n";
	}
}
print MM "i\tj\tvalue\tprb\n";
for (my $i = 0; $i<@MM_X; $i++) {
        print MM $MM_X[$i], "\t", $MM_Y[$i], "\t", $MM_C[$i], "\t", $MM_PRB[$i],"\n";
}
print INSERTIONS "i\tj\tvalue\tprb\n";
for (my $i = 0; $i<=$ins_max_x; $i++) {
	for (my $j = 0; $j<=$ins_max_y; $j++) {
		if (defined($INS[$i][$j])) {
        		print INSERTIONS  $i, "\t", $j, "\t", $INS[$i][$j], "\t", $INS_PRB[$i][$j],"\n";
		}
	}
}

### SET BACKBONE SEQUENCE
my $backbone_length = $end - $begin + 1;
my $count_backbone_chars = 0;
if ($nobackbone == 0) {
	my $q2 = "SELECT base FROM seq_ref where chromosome = $chr and position between $begin and $end order by position";
	my $sth2 = $dbh->prepare($q2); $sth2->execute();
	while (my $res = $sth2->fetchrow_hashref()) {
        	print SEQ $res->{base}, "\n";
		$count_backbone_chars++;
	}
}
for (my $i = $count_backbone_chars; $i < $backbone_length; $i++) {
	print SEQ "-", "\n";
}


close MAPPING;
close INSERTIONS;
close SEQ;

# Plot using R
print STDERR "Auf zu R\n";
#print STDERR "R --slave --vanilla --args '$dir/mapping_4_plotting.txt' $begin $end '$dir/seq_ref.txt' '$dir/mapping_mm.txt' '$dir/mapping_ins.txt' $format < /ebio/abt6/korbinian/pgsp/Plot/Mapping/plot_mapping.R\n";
system("R --slave --vanilla --args '$dir/tmpmap_mapping_4_plotting.txt' $begin $end '$dir/tmpmap_seq_ref.txt' '$dir/tmpmap_mapping_mm.txt' '$dir/tmpmap_mapping_ins.txt' $format $chr < ~/pgsp/Plot/Mapping/plot_mapping.R");
system("rm $dir/tmpmap_mapping_4_plotting.txt $dir/tmpmap_seq_ref.txt $dir/tmpmap_mapping_mm.txt $dir/tmpmap_mapping_ins.txt");

#######################################################

sub parse_seq {
	my ($seq, $x, $y, $prb, $strand) = @_;
	my $len = 0;
	my $nucs = 0;
	my @prb_array;
	my $prb_init = 0;
	for (my $i=0; $i<length($seq); $i++) {
		if (substr($seq, $i, 1) ne "[") {
			$len++;
		
		} else {
			if (substr($seq, $i+1, 1) eq "-") { # Insertion detected
				if ($prb_init == 0) {
					@prb_array = split "", $prb;
					if ($strand eq "P") {
						my @s = reverse @prb_array;
						@prb_array = @s;
					}
					$prb_init = 1;
				}
				$INS[$x+$len-1][$y] .= substr($seq, $i+2, 1);
				$ins_max_x = $x+$len-1 if $x+$len-1 > $ins_max_x;
				$ins_max_y = $y if $y > $ins_max_y;
				$INS_PRB[$x+$len-1][$y] = ord($prb_array[$nucs]) - 50;
			} else {
				if ($prb_init == 0) {
					@prb_array = split "", $prb;
                                        if ($strand eq "P") {
                                                my @s = reverse @prb_array;
                                                @prb_array = @s;
                                        }
					$prb_init = 1;
                                }
				$len++;
				push @MM_C, substr($seq, $i+2, 1);
				push @MM_X, $x+$len-1;
				push @MM_Y, $y;
				push @MM_PRB, ord($prb_array[$nucs]) - 64;
			}
			$i+=3;
		}
		$nucs++;
	}
	return($len);
}

#debug
sub print_mapping_pos {
	my ($x, $last_x) = @_;

	print STDERR "#################################\n";
	print STDERR "x:", $x, "\tlastx:", $last_x, "\n";

	for (my $i = 0; $i < @MAPPING_POS; $i++) {
		print STDERR $i, "\t", $MAPPING_POS[$i], "\n";
	}

}


sub GetCom {
  my @usage = ("\nUsage: $0 --file=<file> --idx=<file> --chr=<int> --begin=<int> --end=<int>

required:
--file\t\tmap.list
--chr\t\tMouse on Mars

One of:
--plotpos\tSingle position with 
or
--begin\t\tbegin position
--end\t\tend position

optional:
--offset\tbp around it (only with plotpos)
--idx\t\tindex of map.list (default: map.list filename concat with .idx)
--database\t\tdatabase (default: solexa)
--format\t\toutput file format (png | pdf) (default: pdf)
-nobackbone\tNo backbone available
\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "file=s", "chr=s", "begin=s", "end=s", "idx=s", "database=s", "format=s", "nobackbone", "plotpos=s", "offset=s");

	die("Please specify an input file\n") unless $CMD{file};
	die("Please specify chr\n") unless $CMD{chr};

	$file = $CMD{file};
	$chr = $CMD{chr};	

	if (defined($CMD{offset})) {
		$offset = $CMD{offset};
        }
        else {
                $offset = 20;
        }

	#############

	$begin = 0;
	$end = 0;

	if (defined($CMD{plotpos})) {
		$begin = $CMD{plotpos} - $offset;
	        $end = $CMD{plotpos} + $offset;
	}
	elsif (defined($CMD{begin})) {
                die("Please specify end\n") unless $CMD{end};
                $begin = $CMD{begin};
                $end = $CMD{end};
        }

	if ($begin == 0 || $end == 0)  {
		die("Please specify (begin, end) or (plotpos, offset)\n") ;
	}

	##############

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

sub connect_to_db {
        my $databaseName = $database;
        my $driver = "mysql";
        my $host = "orb.eb.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}


exit(0);
