#!/usr/bin/perl
#############################################################################
# Author 	Stephan Ossowski, Korbinian Schneeberger
# Date 		05/10/07
# Version	0.95.3.1b
# Function:	Error analysis of reference sequencings 
#############################################################################

use strict;
use warnings;
use Cwd;

### Parse input
my $file_frags           = shift;
my $read_length      = shift;

if( (not defined $file_frags) ) {
	print "Usage example: perl $0 map.list"; 
}

open (FRAGS, $file_frags) or die "Cannot open input file $file_frags\n";

# Phasing and Prephasing
my %lasterror_lastref_vs_readbase = ();
my %lasterror_lastread_vs_readbase = ();
my %lasterror_lastref_vs_refbase = ();
my %lasterror_lastread_vs_refbase = ();
my %error_lastread_vs_readbase = ();
my %error_lastread_vs_refbase = ();
my %error_lastref_vs_readbase = ();
my %error_lastref_vs_refbase = ();
# Base comparison
my %base_comparison = ();
# Errors per position
my %position_error =();
my %position_occ = ();
# Errors per read
my %read_error_number = ();
# Cumulation of errors
my @error_dist_2 = ();
my @error_dist_3 = ();
my @error_dist_4 = ();
my @error_dist_max = ();
# Occurrence of first error in read
my @error_first_occurrence = ();

my $count_reads = 0;

### Scan input file
while( <FRAGS> ) {
	chomp($_);
	my @members = split(/\t/, $_);
	my $chr            = $members[0];
	my $position       = $members[1];
	my $s              = $members[2];
	my $solexa_id      = $members[3];
	my $orientation    = $members[4];
	my $mismatch       = $members[5];
	my $hits           = $members[6];
	my $used_length    = $members[7];
	my $offset         = $members[8];
	my $prb_line       = $members[9];
	my $qval_line	   = $members[10];
	my $chas_line 	   = $members[11];

	$count_reads++;
	print STDERR $count_reads, "\n" if $count_reads % 1000000 == 0;

	if ($hits == 1) { 

		### Rev. complement read if mapped palindromic ###################################################################
		if ($orientation eq "P") {
			my $t = reverse $s;
			$t =~ s/\]/\+/g;
			$t =~ s/\[/\]/g;
			$t =~ s/\+/\[/g;
			$t =~ tr/ACGTacgt/TGCAtgca/;
			# Swap bases between brackets
			for (my $i=0; $i<length($t); $i++) {
				if (substr($t, $i, 1) eq "[") {
					my $r = substr($t, 0, $i+1).substr($t, $i+2, 1).substr($t, $i+1, 1).substr($t, $i+3, length($t)-($i+3));
					$t = $r;
					$i += 2;
				}
			}
			$s = $t;
		}

	
		my $pos_in_read = 0;
		my $first_occ = -1;
		my $last_readbase;
		my $last_refbase;
		my $last_errorflag = 0;
		my $errors_per_read = 0;
		my $i = 0;
		my @error_occ = ();

	        while($i < length($s)) {

			my $readbase;
			my $refbase;  
			my $errorflag;              
	
			### Scan current read position ####
                	if(substr($s, $i, 1) ne "[") {
				$readbase = substr($s, $i, 1);
				$refbase = $readbase;
				$pos_in_read++;
				$errorflag = 0;
				$i++;
			}
			else {
                	        $i++; $refbase = substr($s, $i, 1);
				$i++; $readbase = substr($s, $i, 1);

        	                if ($readbase ne "-") {
					$pos_in_read++;
                        	}

				$errorflag = 1;
				$i += 2;
	                }

			### Base comparison ###
			$base_comparison{$refbase.$readbase}++;

			### Phasing and Prephasing ###
			if ($last_errorflag == 1) {
				if (defined($last_refbase) and $readbase ne "-") {
					$lasterror_lastref_vs_readbase{$last_refbase.$readbase}++;
				}
				if (defined($last_readbase) and $readbase ne "-") {
	        	                $lasterror_lastread_vs_readbase{$last_readbase.$readbase}++;
				}
				if (defined($last_refbase) and $refbase ne "-") {
        	                        $lasterror_lastref_vs_refbase{$last_refbase.$refbase}++;
                	        }
                        	if (defined($last_readbase) and $refbase ne "-") {
	                                $lasterror_lastread_vs_refbase{$last_readbase.$refbase}++;
        	                }
			}

			if ($errorflag == 1) {
				if (defined($last_readbase) and $readbase ne "-") {
					$error_lastread_vs_readbase{$last_readbase.$readbase}++ 
				}
				if (defined($last_readbase) and $refbase ne "-") {
					$error_lastread_vs_refbase{$last_readbase.$refbase}++ 	
				}
				if (defined($last_refbase) and $readbase ne "-") {
        	                        $error_lastref_vs_readbase{$last_refbase.$readbase}++
                	        }
                        	if (defined($last_refbase) and $refbase ne "-") {
                                	$error_lastref_vs_refbase{$last_refbase.$refbase}++
	                        }
			}

			### Error per read position
			if ($errorflag == 1) {
				$position_error{$pos_in_read}++;
				push @error_occ, $pos_in_read;
			}
			$position_occ{$pos_in_read}++;

			### Errors per read
			if ($errorflag == 1) {
				$errors_per_read++;
			}

			### First occurrence of an error
			if ($errorflag == 1 and $first_occ == -1) {
				push @error_first_occurrence, $pos_in_read;	
				$first_occ = $pos_in_read;
			}

			### Reset variables ###
			$last_readbase = $readbase if $readbase ne "-";
                	$last_refbase = $refbase if $refbase ne "-";
			$last_errorflag = $errorflag;
        	}

		$read_error_number{$errors_per_read}++;
		
		# Error distances
		if ($errors_per_read == 2) {
			push @error_dist_2, ($error_occ[1] - $error_occ[0]);
			push @error_dist_max, ($error_occ[1] - $error_occ[0]);	
		}
		if ($errors_per_read == 3) {
			push @error_dist_3, ($error_occ[1] - $error_occ[0]);
			push @error_dist_3, ($error_occ[2] - $error_occ[0]);
			push @error_dist_3, ($error_occ[2] - $error_occ[1]);
			push @error_dist_max, ($error_occ[2] - $error_occ[0]);
		}
		if ($errors_per_read == 4) {
			push @error_dist_4, ($error_occ[1] - $error_occ[0]);
			push @error_dist_4, ($error_occ[2] - $error_occ[0]);
			push @error_dist_4, ($error_occ[3] - $error_occ[0]);
			push @error_dist_4, ($error_occ[2] - $error_occ[1]);
			push @error_dist_4, ($error_occ[3] - $error_occ[1]);
			push @error_dist_4, ($error_occ[3] - $error_occ[2]);
			push @error_dist_max, ($error_occ[3] - $error_occ[0]);
		}


	}

}

### Create random distance distributions ###
my @error_dist_max_rand = ();

my @error_dist_2_rand = ();
for (my $i = 0; $i<@error_dist_2; $i++) {
	my @err = ();
	$err[0] = int(rand(35)+1);
	$err[1] = $err[0];
	while ($err[0] == $err[1]) {
		$err[1] = int(rand(35)+1);
	}
	my @err_sort = sort {$a <=> $b} @err;
	push @error_dist_2_rand, $err_sort[1] - $err_sort[0];
	push @error_dist_max_rand, $err_sort[1] - $err_sort[0];
}

my @error_dist_3_rand = ();
for (my $i = 0; $i<@error_dist_3/3; $i++) { # 3 distance per loop
	my @err = ();
        $err[0] = int(rand(35)+1);
        $err[1] = $err[0];
        while ($err[0] == $err[1]) {
                $err[1] = int(rand(35)+1);
        }
	$err[2] = $err[1];
	while ($err[2] ==  $err[1] or $err[2] == $err[0]) {
		$err[2] = int(rand(35)+1);
	}
        my @err_sort = sort {$a <=> $b} @err;
        push @error_dist_3_rand, $err_sort[1] - $err_sort[0];
	push @error_dist_3_rand, $err_sort[2] - $err_sort[0];
	push @error_dist_3_rand, $err_sort[2] - $err_sort[1];
	push @error_dist_max_rand, $err_sort[2] - $err_sort[0];
}	

my @error_dist_4_rand = ();
for (my $i = 0; $i<@error_dist_4/6; $i++) { # 6 mm per loop
	my @err = ();
	$err[0] = int(rand(35)+1);
        $err[1] = $err[0];
        while ($err[0] == $err[1]) {
                $err[1] = int(rand(35)+1);
        }
        $err[2] = $err[1];
        while ($err[2] ==  $err[1] or $err[2] == $err[0]) {
                $err[2] = int(rand(35)+1);
        }
	$err[3] = $err[2];
	while ($err[3] ==  $err[1] or $err[3] == $err[0] or $err[3] == $err[2]) {
		$err[3] = int(rand(35)+1);
	}
	my @err_sort = sort {$a <=> $b} @err;
	push @error_dist_4_rand, $err_sort[1] - $err_sort[0];
        push @error_dist_4_rand, $err_sort[2] - $err_sort[0];
        push @error_dist_4_rand, $err_sort[2] - $err_sort[1];
	push @error_dist_4_rand, $err_sort[3] - $err_sort[2];
        push @error_dist_4_rand, $err_sort[3] - $err_sort[0];
        push @error_dist_4_rand, $err_sort[3] - $err_sort[1];
	push @error_dist_max_rand, $err_sort[3] - $err_sort[0];
}

### Print prephasing
# ERRORs compared with the upfollowing base
open FILE, "> error_prephasing_ref_vs_read.txt";
foreach my $key (sort {$a cmp $b} keys %lasterror_lastref_vs_readbase) { 
	print FILE substr($key, 0, 1), "\t", substr($key, 1, 1), "\t", $lasterror_lastref_vs_readbase{$key}, "\n";
}
close FILE;
open FILE, "> error_prephasing_read_vs_read.txt";
foreach my $key (sort {$a cmp $b} keys %lasterror_lastread_vs_readbase) {
	print FILE substr($key, 0, 1), "\t", substr($key, 1, 1), "\t", $lasterror_lastread_vs_readbase{$key}, "\n";
}
close FILE;
open FILE, "> error_prephasing_ref_vs_ref.txt";
foreach my $key (sort {$a cmp $b} keys %lasterror_lastref_vs_refbase) {
	print FILE substr($key, 0, 1), "\t", substr($key, 1, 1), "\t", $lasterror_lastref_vs_refbase{$key}, "\n";
}
close FILE;
open FILE, "> error_prephasing_read_vs_ref.txt";
foreach my $key (sort {$a cmp $b} keys %lasterror_lastread_vs_refbase) {
	print FILE substr($key, 0, 1), "\t", substr($key, 1, 1), "\t", $lasterror_lastread_vs_refbase{$key}, "\n";
}
close FILE;

### Print phasing
# ERRORs compared with the previous base:
open FILE, "> error_phasing_ref_vs_read.txt";
foreach my $key (sort {$a cmp $b} keys %error_lastref_vs_readbase) {
        print FILE substr($key, 0, 1), "\t", substr($key, 1, 1), "\t", $error_lastref_vs_readbase{$key}, "\n";
}
close FILE;
open FILE, "> error_phasing_read_vs_read.txt";
foreach my $key (sort {$a cmp $b} keys %error_lastread_vs_readbase) {
        print FILE substr($key, 0, 1), "\t", substr($key, 1, 1), "\t", $error_lastread_vs_readbase{$key}, "\n";
}
close FILE;
open FILE, "> error_phasing_ref_vs_ref.txt";
foreach my $key (sort {$a cmp $b} keys %error_lastref_vs_refbase) {
        print FILE substr($key, 0, 1), "\t", substr($key, 1, 1), "\t", $error_lastref_vs_refbase{$key}, "\n";
}
close FILE;
open FILE, "> error_phasing_read_vs_ref.txt";
foreach my $key (sort {$a cmp $b} keys %error_lastread_vs_refbase) {
        print FILE substr($key, 0, 1), "\t", substr($key, 1, 1), "\t", $error_lastread_vs_refbase{$key}, "\n";
}
close FILE;

### Print Base comparison
open FILE, "> error_basecomparison.txt";
foreach my $key (sort {$a cmp $b} keys %base_comparison) {
	print FILE substr($key, 0, 1), "\t", substr($key, 1, 1), "\t", $base_comparison{$key}, "\n";
}
close FILE;

### Print Position Errors
open FILE, "> error_positionwise.txt";
foreach my $key (sort {$a <=> $b} keys %position_error) {
	print FILE $key, "\t", $position_error{$key}, "\t", $position_occ{$key}, "\t", $position_error{$key}/$position_occ{$key}, "\n";
}
close FILE;

### Print error number per read
open FILE, "> error_num_per_read.txt";
foreach my $key (sort {$a <=> $b} keys %read_error_number) {
	print FILE $key, "\t", $read_error_number{$key}, "\n";
}
close FILE;

### Print first error occurrence distribution
open FILE, "> error_first_occurrence.txt";
foreach my $key (@error_first_occurrence) {
        print FILE $key, "\n";
}
close FILE;

### Print error distances 
open FILE, "> error_distances_2.txt";
foreach my $dist (@error_dist_2) {
	print FILE $dist, "\n";
}
close FILE;
open FILE, "> error_distances_3.txt";
foreach my $dist (@error_dist_3) {
	print FILE $dist, "\n";
}
close FILE;
open FILE, "> error_distances_4.txt";
foreach my $dist (@error_dist_4) {
        print FILE $dist, "\n";
}
close FILE;
open FILE, "> error_distances_max.txt";
foreach my $dist (@error_dist_max) {
	print FILE $dist, "\n";
}
close FILE;

### Print random error distances
open FILE, "> error_distances_2_rand.txt";
foreach my $dist (@error_dist_2_rand) {
        print FILE $dist, "\n";
}
close FILE;
open FILE, "> error_distances_3_rand.txt";
foreach my $dist (@error_dist_3_rand) {
        print FILE $dist, "\n";
}
close FILE;
open FILE, "> error_distances_4_rand.txt";
foreach my $dist (@error_dist_4_rand) {
        print FILE $dist, "\n";
}
close FILE;
open FILE, "> error_distances_max_rand.txt";
foreach my $dist (@error_dist_max_rand) {
        print FILE $dist, "\n";
}
close FILE;

exit(0);



