#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Cwd;

my $dir;
my %CMD;
my $min;
my $max;
my $mid;

my $PAIR_UNPAIRED_SINGLE_SINGLE;
my $PAIR_UNPAIRED_SINGLE_MULTIPLE;
my $PAIR_UNPAIRED_MULTIPLE_SINGLE;
my $PAIR_UNPAIRED_MULTIPLE_MULTIPLE;
my $PAIR_PAIRED_SINGLE_SINGLE;
my $PAIR_PAIRED_MULTIPLE_MULTIPLE;
my $SINGLETON_SINGLE_1;
my $SINGLETON_MULTIPLE_1;
my $SINGLETON_SINGLE_2;
my $SINGLETON_MULTIPLE_2;

GetCom();

#### Log to commaend line options

open LOG, "> mappingcorrection.txt";
foreach my $key (keys %CMD) {
	print LOG $key, "\t", $CMD{$key}, "\n";
}
close LOG;

######## Collect map.lists
my $cwd = cwd();

chdir($cwd."/".$dir."/1");
my $file_string = "cat ";
foreach my $length (glob("length_*")) {
	$file_string .= $length."/map.list.sorted_id ";
}
system($file_string." > map.list.sorted_id");

chdir($cwd."/".$dir."/2");
$file_string = "cat ";
foreach my $length (glob("length_*")) {
        $file_string .= $length."/map.list.sorted_id ";
}
system($file_string." > map.list.sorted_id");

chdir($cwd."/".$dir);
###########################

open MAP1, "1/map.list.sorted_id" or die "Cannot open inputfile\n";
open MAP2, "2/map.list.sorted_id" or die "Cannot open inputfile\n";

open MAP1CORR, "> 1/map.list.1.corrected.sorted_id" or die "Cannot open outputfile";
open MAP2CORR, "> 2/map.list.2.corrected.sorted_id" or die "Cannot open outputfile";
open MAP_JOIN, "> pair.list.unsorted" or die "Cannot open outputfile";
open MAP_DIST, "> insert_size.txt" or die "Cannot open outputfile\n";

open MAP_STAT, "> mappingcorrection.stat";

my $id1 = "";
my $id2 = "";

my $id1_linebuffer = "";
my $id2_linebuffer = "";

my @id1_lines = ();
my @id2_lines = ();

my %READS1 = ();
my %READS2 = ();

my $count = 0;
while (not (eof(MAP1) and eof(MAP2))) {

	$count++;
	#print STDERR $count, "\t>", $id1, "<\t>", $id2, "<\n";
	print $count, "\n" if $count%10000 == 0;	


	###################################################################
	### Read in the next id from both files  ##########################
	###################################################################


	###################################################################
	# Read in lines of read in map.list.unsorted.1 if not already done

	if ($id1 eq "") {
		if ($id1_linebuffer ne "") {
			my @entries = split " ", $id1_linebuffer;
			$id1 = $entries[3];
        		push @id1_lines, $id1_linebuffer;
			$id1_linebuffer = "";
	        }

		FILE: while (my $line = <MAP1>) {
			my @entries = split " ", $line;
			my $id = $entries[3];

			if ($id1 eq "") {
				$id1 = $id;
				push @id1_lines, $line;
			}
			else {
				if ($id1 == $id) {
					push @id1_lines, $line;
				}
				else {
					$id1_linebuffer = $line;
					last FILE;
				}
			}
		}
	}

	###################################################################
	# Read in lines of read in map.list.unsorted.2 if not already done

        if ($id2 eq "") {
                if ($id2_linebuffer ne "") {
                        my @entries = split " ", $id2_linebuffer;
                        $id2 = $entries[3];
                        push @id2_lines, $id2_linebuffer;
                        $id2_linebuffer = "";
                }

                FILE: while (my $line = <MAP2>) {
                        my @entries = split " ", $line;
                        my $id = $entries[3];

                        if ($id2 eq "") {
                                $id2 = $id;
                                push @id2_lines, $line;
                        }
                        else {
                                if ($id2 == $id) {
                                        push @id2_lines, $line;
                                }
                                else {
                                        $id2_linebuffer = $line;
                                        last FILE;
                                }
                        }
                }
        }

	################################################################
	#### Compare mappings of paired reads if not paired flush ######
	################################################################
	#### Both reads were mapped                               ######
	################################################################
	if ($id1 == $id2) {

		####################################
		### Check for paired mappings

		my @pair1_partner = ();
		my @pair2_partner = ();
		my @pair_distance = ();

		my %pairedhits1_count = ();
		my %pairedhits2_count = ();

		my $pairings = 0;

		for (my $i = 0; $i<@id1_lines; $i++) {
			my @entries1 = split " ", $id1_lines[$i];
			my $chr1 = $entries1[0];
			my $pos1 = $entries1[1];
			my $dir1 = $entries1[4];
			my $len1 = $entries1[7];

			for (my $j = 0; $j<@id2_lines; $j++) {
				my @entries2 = split " ", $id2_lines[$j];
				my $chr2 = $entries2[0];
	                        my $pos2 = $entries2[1];
                	        my $dir2 = $entries2[4];
				my $len2 = $entries2[7];
				my $paired = 0;
		
				# Are the mappings paired?
				my $dist;
				if ($dir1 ne $dir2 and $chr1 == $chr2) {
					if ($dir1 eq "D") {
						$dist = ($pos2+$len2)-$pos1;
					}
					else {
						$dist = ($pos1+$len1)-$pos2;
					}
				
					if ($dist >= $min and $dist <= $max) {	
						$paired = 1;
					}					
				}

				# Store paired mappings
				if ($paired == 1) {

					$pair1_partner[$pairings] = $i;
					$pair2_partner[$pairings] = $j;
					$pair_distance[$pairings] = $dist;

					$pairedhits1_count{$i}++;
					$pairedhits2_count{$j}++;

					$pairings++;
				}
			}
		}	

		#############################################
		### No pair between mapping could be found
		if ($pairings == 0) {

			# Set Pairing flag
			my $pair_flag = "A";

			# Set repetitive rep_flags
			my $rep_flag1;
                        my $rep_flag2;
			if (@id1_lines == 1 and @id2_lines == 1) {
                                $PAIR_UNPAIRED_SINGLE_SINGLE++;
				$rep_flag1 = "U"; $rep_flag2 = "U";
                        }
                        elsif (@id1_lines == 1 and @id2_lines > 1) {
                                $PAIR_UNPAIRED_SINGLE_MULTIPLE++;
				$rep_flag1 = "U"; $rep_flag2 = "R";
                        }
                        elsif (@id1_lines > 1 and @id2_lines == 1) {
                                $PAIR_UNPAIRED_MULTIPLE_SINGLE++;
				$rep_flag1 = "R"; $rep_flag2 = "U";
                        }
                        elsif (@id1_lines > 1 and @id2_lines > 1) {
                                $PAIR_UNPAIRED_MULTIPLE_MULTIPLE++;
				$rep_flag1 = "R"; $rep_flag2 = "R";
                        }

			# Print 
			for (my $i = 0; $i<@id1_lines; $i++) {
				chomp($id1_lines[$i]);
                       	        print MAP1CORR $id1_lines[$i], "\t1\n";

				my @tmp1 = split " ", $id1_lines[$i];
				my $chr1 = $tmp1[0];
				my $pos1 = $tmp1[1];
				my $id = $tmp1[3];
				my $dir1 = $tmp1[4];
				my $len1 = $tmp1[7];

				print MAP_JOIN $chr1, "\t", $pos1, "\t", $len1, "\t", $dir1, "\t1\t", $id, "\t", @id1_lines+0, "\t", $pair_flag, $rep_flag1, $rep_flag2;
				for (my $j=0; $j<@id2_lines; $j++) {
					my @tmp2 = split " ", $id2_lines[$j];
					my $chr2 = $tmp2[0];
	                                my $pos2 = $tmp2[1];
                	                my $dir2 = $tmp2[4];
					my $len2 = $tmp2[7];
					print MAP_JOIN "\t", $chr2, "\t", $pos2, "\t", $len2, "\t", $dir2;
				}
				print MAP_JOIN "\t", @id2_lines+0, "\n";
                        }
			
			for (my $i = 0; $i<@id2_lines; $i++) {
				chomp($id2_lines[$i]);
                                print MAP2CORR $id2_lines[$i], "\t2\n";

                                my @tmp1 = split " ", $id2_lines[$i];
                                my $chr1 = $tmp1[0];
                                my $pos1 = $tmp1[1];
                                my $id = $tmp1[3];
                                my $dir1 = $tmp1[4];
				my $len1 = $tmp1[7];

                                print MAP_JOIN $chr1, "\t", $pos1, "\t", $len1, "\t", $dir1, "\t2\t", $id, "\t", @id1_lines+0, "\t", $pair_flag, $rep_flag2, $rep_flag1;
                                for (my $j=0; $j<@id1_lines; $j++) {
                                        my @tmp2 = split " ", $id1_lines[$j];
                                        my $chr2 = $tmp2[0];
                                        my $pos2 = $tmp2[1];
                                        my $dir2 = $tmp2[4];
					my $len2 = $tmp2[7];
                                        print MAP_JOIN "\t", $chr2, "\t", $pos2, "\t", $len2, "\t", $dir2;
                                }
                                print MAP_JOIN "\t", @id1_lines+0, "\n";
                        }
			



		}
                ######################################################
                ### Pairs could be found
		else {
			######################################################################################
			# Resolve double pairs -- optimal solution too expensive for this artificail problem

			#print "#############\n";
			#for (my $pair=0; $pair < @pair1_partner; $pair++) {
			#	print $pair1_partner[$pair], "\t", $pair2_partner[$pair], "\n";
			#}
			#print "+++++++++++++\n";

			for (my $pair=0; $pair < @pair1_partner; $pair++) { ## speed up would be to compare the nb of different entries to length of @pair1_partner
				my $i = $pair1_partner[$pair];
				my $j = $pair2_partner[$pair];
				my $dist = $pair_distance[$pair];
				
				if (not($pairedhits1_count{$i} == 1 and $pairedhits2_count{$j} == 1)) {
					# run through all other hits and delete one if using same mappings 
					for (my $k = $pair+1; $k<@pair1_partner; $k++) {
						if ($i == $pair1_partner[$k] or $j == $pair2_partner[$k]) {
							if (abs($dist-$mid) <= abs($pair_distance[$k]-$mid)) {
								$pairedhits1_count{$pair1_partner[$k]}--;
                                                                $pairedhits2_count{$pair2_partner[$k]}--;

								splice(@pair1_partner, $k, 1);
								splice(@pair2_partner, $k, 1);
								splice(@pair_distance, $k, 1);
			
								$k--;
							}
							else {
								$pairedhits1_count{$i}--;
                                                                $pairedhits2_count{$j}--;

								splice(@pair1_partner, $pair, 1);
                                                                splice(@pair2_partner, $pair, 1);
                                                                splice(@pair_distance, $pair, 1);

								$pair--;
								$k = @pair1_partner+0; 
							}
						}
					}	
				}
			}

			#######################################################################################
			### Print out 

			my $pair_flag = "P";

			my $rep_flag1;
                        my $rep_flag2;
                        if (@pair1_partner == 1 and @pair2_partner == 1) {
                                $PAIR_PAIRED_SINGLE_SINGLE++;
                                $rep_flag1 = "U"; $rep_flag2 = "U";

				print MAP_DIST $pair_distance[0], "\n";
                        }
                        else  {
                                $PAIR_PAIRED_MULTIPLE_MULTIPLE++;
				$rep_flag1 = "R"; $rep_flag2 = "R";
			}

			for (my $pair = 0; $pair < @pair1_partner; $pair++) {
				my @tmp1 = split " ", $id1_lines[$pair1_partner[$pair]];
				my @tmp2 = split " ", $id2_lines[$pair2_partner[$pair]];
				my $id = $tmp1[3];
				my $chr = $tmp1[0];
				my $hits = @pair1_partner+0;

				my $pos1 = $tmp1[1];
				my $dir1 = $tmp1[4];
				my $len1 = $tmp1[7];

				my $pos2 = $tmp2[1];
				my $dir2 = $tmp2[4];
				my $len2 = $tmp2[7];

				# debug 
				if ($tmp1[0] != $tmp2[0] or @pair1_partner+0 != @pair2_partner+0) {
					print STDERR "Pair is not uniform\n";
					exit(1);
				}

				print MAP1CORR $tmp1[0], "\t", $tmp1[1], "\t", $tmp1[2], "\t", $tmp1[3], "\t", $tmp1[4], "\t", $tmp1[5], "\t", @pair1_partner+0, "\t";
				print MAP1CORR $tmp1[7], "\t", $tmp1[8], "\t", $tmp1[9], "\t", $tmp1[10], "\t", $tmp1[11], "\t", "1", "\n";
				print MAP2CORR $tmp2[0], "\t", $tmp2[1], "\t", $tmp2[2], "\t", $tmp2[3], "\t", $tmp2[4], "\t", $tmp2[5], "\t", @pair2_partner+0, "\t";
                                print MAP2CORR $tmp2[7], "\t", $tmp2[8], "\t", $tmp2[9], "\t", $tmp2[10], "\t", $tmp2[11], "\t", "2", "\n";
				
				print MAP_JOIN $chr, "\t", $pos1, "\t", $len1, "\t", $dir1, "\t1\t", $id, "\t", $hits, "\t", $pair_flag, $rep_flag2, $rep_flag1, "\t", $chr, "\t", $pos2, "\t", $len2, "\t", $dir2, "\t", $hits, "\n";
				print MAP_JOIN $chr, "\t", $pos2, "\t", $len2, "\t", $dir2, "\t2\t", $id, "\t", $hits, "\t", $pair_flag, $rep_flag1, $rep_flag2, "\t", $chr, "\t", $pos1, "\t", $len1, "\t", $dir1, "\t", $hits, "\n";
			}

		}

		$id1 = "";
                @id1_lines = ();
		$id2 = "";
                @id2_lines = ();

	}
	####################################################
	#### Only one read was mapped ######################
	####################################################
	else {
		if ($id2 gt $id1) {
			###############
			# flush id1
	
			# Set pairing flag
			my $pair_flag = "M";

			# Set repetitive flag 
			my $rep_flag1;
			my $rep_flag2 = "M";

			if (@id1_lines == 1) {
				$SINGLETON_SINGLE_1++;
				$rep_flag1 = "U";
			}
			else {
				$SINGLETON_MULTIPLE_1++;
				$rep_flag1 = "R";
			}

			# Print 
			for (my $i = 0; $i<@id1_lines; $i++) {
				chomp($id1_lines[$i]);
				print MAP1CORR $id1_lines[$i], "\t", "1", "\n";
				
				my @tmp = split " ", $id1_lines[$i];
				my $chr = $tmp[0];
				my $pos = $tmp[1];
				my $id = $tmp[3];
				my $dir = $tmp[4];
				my $len = $tmp[7];
				print MAP_JOIN $chr, "\t", $pos, "\t", $len, "\t", $dir, "\t1\t", $id, "\t", @id1_lines+0, "\t", $pair_flag, $rep_flag1, $rep_flag2, "\n";
			}
			
			$id1 = "";
                        @id1_lines = ();	
		}
		else {
			###############
			# flush id2

			# Set pairing flag
			my $pair_flag = "M";

			# Set repetitive flag
			my $rep_flag1;
                        my $rep_flag2 = "M";

                        if (@id2_lines == 1) {
                                $SINGLETON_SINGLE_2++;
                                $rep_flag1 = "U";
                        }
                        else {
                                $SINGLETON_MULTIPLE_2++;
                                $rep_flag1 = "R";
                        }

			# Print
			for (my $i = 0; $i<@id2_lines; $i++) {
				chomp($id2_lines[$i]);
                                print MAP2CORR $id2_lines[$i], "\t", "2", "\n";

                                my @tmp = split " ", $id2_lines[$i];
                                my $chr = $tmp[0];
                                my $pos = $tmp[1];
                                my $id = $tmp[3];
                                my $dir = $tmp[4];
				my $len = $tmp[7];
                                print MAP_JOIN $chr, "\t", $pos, "\t", $len, "\t", $dir, "\t2\t", $id, "\t", @id2_lines+0, "\t", $pair_flag, $rep_flag1, $rep_flag2, "\n";
                        }

                        $id2 = "";
                        @id2_lines = ();

		}

	}

}

close MAP1CORR;
close MAP2CORR;
close MAP_JOIN;


#### Sort by genomic position
system("sort -n -k1 -k2 1/map.list.1.corrected.sorted_id > 1/map.list.1");
system("sort -n -k1 -k2 2/map.list.2.corrected.sorted_id > 2/map.list.2");
system("sort -n -k1 -k2 pair.list.unsorted > pair.list");

system("rm pair.list.unsorted 1/map.list.sorted_id 2/map.list.sorted_id 1/map.list.1.corrected.sorted_id 2/map.list.2.corrected.sorted_id");

### Statistics

close MAP_DIST;
system("R --slave --vanilla --args '$dir/insert_size.txt' $min $max < /ebio/abt6/korbinian/pgsp/Plot/PE/plot_insert_dist.R");


print MAP_STAT "*****************************************************************\n";
print MAP_STAT "Missing partner mapping:\t", ($SINGLETON_SINGLE_1+$SINGLETON_MULTIPLE_1+$SINGLETON_SINGLE_2+$SINGLETON_MULTIPLE_2), "\n";
print MAP_STAT "\tSingletons:\t\t", ($SINGLETON_SINGLE_1+$SINGLETON_SINGLE_2), " (", $SINGLETON_SINGLE_1, " + ", $SINGLETON_SINGLE_2, ")\n";
print MAP_STAT "\tMultiple:\t\t", ($SINGLETON_MULTIPLE_1+$SINGLETON_MULTIPLE_2), " (",$SINGLETON_MULTIPLE_1, " + ",  $SINGLETON_MULTIPLE_2, ")\n";

print MAP_STAT "\n*****************************************************************\n";
print MAP_STAT "Pairs with disturbed mappings:\t", ($PAIR_UNPAIRED_SINGLE_SINGLE+$PAIR_UNPAIRED_SINGLE_MULTIPLE+$PAIR_UNPAIRED_MULTIPLE_SINGLE+$PAIR_UNPAIRED_MULTIPLE_MULTIPLE), "\n";
print MAP_STAT "\tSS:\t".$PAIR_UNPAIRED_SINGLE_SINGLE."\n";
print MAP_STAT "\tMM:\t".$PAIR_UNPAIRED_MULTIPLE_MULTIPLE."\n";
print MAP_STAT "\tSM:\t".$PAIR_UNPAIRED_SINGLE_MULTIPLE."\n";
print MAP_STAT "\tMS:\t".$PAIR_UNPAIRED_MULTIPLE_SINGLE."\n";

print MAP_STAT "\n*****************************************************************\n";
print MAP_STAT "Pairs:\t", ($PAIR_PAIRED_SINGLE_SINGLE+$PAIR_PAIRED_MULTIPLE_MULTIPLE), "\n";
print MAP_STAT "\tSS:\t".$PAIR_PAIRED_SINGLE_SINGLE."\n";
print MAP_STAT "\tMM:\t".$PAIR_PAIRED_MULTIPLE_MULTIPLE."\n";


close MAP_STAT;

exit(0);


sub GetCom {

  my @usage = ("
Usage: 	$0 

--min\t\tMinimal insert size
--max\t\tMaximal insert size
--dir\t\tLane folder

\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, 
			"min=s", "max=s", "dir=s");

	die("Please specify max insert size\n") unless defined($CMD{max});
	die("Please specify minimum insert size\n") unless defined($CMD{min});
	die("Please specify lane folder\n") unless defined($CMD{dir});

        $max = $CMD{max};
	$min = $CMD{min};
	$mid = (($max - $min)/2) + $min;
	$dir = $CMD{dir};

}

