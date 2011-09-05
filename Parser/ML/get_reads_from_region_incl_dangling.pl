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
#  Module: Parser::ML::get_reads_from_region_incl_dangling.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 map.list missing_left_over.fl regions outputFolder
region file = <chr> <begin> <end>
Needs to be sorted!\n";
my $map_list = shift or die $usage;
my $mlo = shift or die $usage;
my $regions = shift or die $usage;
my $folder = shift or die $usage;

if (-e $folder) { die("Folder $folder exists\n"); }
system("mkdir -p $folder");

### Read left over reads
my %MLO = ();
open MLOF, $mlo or die "cannot open $mlo\n";
while(<MLOF>) {
        my @a = split " ", $_;
        $MLO{$a[0]} = $a[1];
}
close MLOF;

### Read target regions
my @REG_CHR = ();
my @REG_START = ();
my @REG_END = ();
open REG, $regions or die "cannot open $regions\n";
while(<REG>) {
	my @a = split " ", $_;
	push @REG_CHR, $a[0];
	push @REG_START, $a[1];
	push @REG_END, $a[2];
}
close REG;

### Get reads from target region
my $in_flag = 0;

my %first_reads = ();
my %second_reads = ();
my %first_missing = ();
my %second_missing = ();
my %singletons = ();

my %first_locations = ();
my %second_locations = ();
my %singleton_locations = ();

my $assem_id = 0;

open LOG, ">$folder/region.log.txt" or die "cannot open log file\n";
open FILE, $map_list or die "cannnot open $map_list\n";
my $fp_last = tell(FILE);

while (@REG_CHR+0 > 0) {

	my $fp_line = tell(FILE);
	my $l = <FILE>;
	my @a = split " ", $l;
	my $chr = $a[0];
	my $pos = $a[1];
	my $rlength = $a[7]; # ???
	my $pe = $a[9]; #??

	if ($in_flag) {
		# still in?
		if ($pos+$rlength <= $REG_END[0]) {
			set_read($l);
		}
		else { # not in anymore


			# write out reads
			foreach my $id (keys %first_reads) {
				if (defined($second_reads{$id})) {
					print PAIR ">", $id, "\n", $first_reads{$id}, "\n";
					print PAIR ">", $id, "\n", $second_reads{$id}, "\n";
					delete $second_reads{$id};
					print LOG $id, "\t1\t", $REG_CHR[0], "\t", $first_locations{$id}, "\n";
					print LOG $id, "\t2\t", $REG_CHR[0], "\t", $second_locations{$id}, "\n";
				}
				else {
					print SINGLETON ">", $id, "\n", $first_reads{$id}, "\n";
					print LOG $id, "\tS\t", $REG_CHR[0], "\t", $first_locations{$id}, "\n";
				}
				delete $first_reads{$id};
			}
			foreach my $id (keys %second_reads) {
				print SINGLETON ">", $id, "\n", $second_reads{$id}, "\n";
				delete $second_reads{$id};
				print LOG $id, "\tS\t", $REG_CHR[0], "\t", $second_locations{$id}, "\n";
			}
			foreach my $id (keys %singletons) {
				print SINGLETON ">", $id, "\n", $singletons{$id}, "\n";
				delete $singletons{$id};
				print LOG $id, "\tS\t", $REG_CHR[0], "\t", $singleton_locations{$id}, "\n";
			}
			foreach my $id (keys %first_missing) {
				if (defined($MLO{$id})) {
					print PAIR ">", $id, "\n", $first_missing{$id}, "\n";
					print PAIR ">", $id, "\n", $MLO{$id}, "\n";
					delete $MLO{$id};
					delete $first_missing{$id};
					print LOG $id, "\t1\t", $REG_CHR[0], "\t", $first_locations{$id}, "\n";
					print LOG $id, "\t2\t", $REG_CHR[0], "\t", "M", "\n";
				}
				else {
					print SINGLETON ">", $id, "\n", $first_missing{$id}, "\n";
					delete $first_missing{$id};
                                        print LOG $id, "\tS\t", $REG_CHR[0], "\t", $first_locations{$id}, "\n";
				}
			}
			foreach my $id (keys %second_missing) {
                                if (defined($MLO{$id})) {
                                        print PAIR ">", $id, "\n", $second_missing{$id}, "\n";
                                        print PAIR ">", $id, "\n", $MLO{$id}, "\n";
                                        delete $MLO{$id};
                                        delete $second_missing{$id};
                                        print LOG $id, "\t1\t", $REG_CHR[0], "\t", "M", "\n";
					print LOG $id, "\t2\t", $REG_CHR[0], "\t", $second_locations{$id}, "\n";
                                }
                                else {
                                        print SINGLETON ">", $id, "\n", $second_missing{$id}, "\n";
                                        delete $second_missing{$id};
					print LOG $id, "\tS\t", $REG_CHR[0], "\t", $second_locations{$id}, "\n";
                                }
                        }
				

			# delete first region
			shift @REG_CHR;
			shift @REG_START;
			shift @REG_END;

			# Yet another region?
			if (@REG_CHR+0 == 0) {
				exit(0);
			}

			# if too far in file => restart
			if ($REG_CHR[0] < $chr or ($REG_CHR[0] == $chr and $REG_START[0] < $pos)) {
				#close FILE;
				#open FILE, $map_list or die "cannnot open $map_list\n";
				seek(FILE, $fp_last, 0); 
			}
		
			$in_flag = 0;
		}
	}
	else {
		# now in?
		if ($a[0] == $REG_CHR[0] and $a[1] >= $REG_START[0]) {
			set_read($l);
			$in_flag = 1;

			# Set filepointer to jump back to
			$fp_last = $fp_line;

			# open files
			close PAIR;
			close SINGLETON;

			$assem_id++;

			print LOG "#$assem_id $REG_CHR[0] $REG_START[0] $REG_END[0]\n";

			open PAIR, ">$folder/assem_id$assem_id.p.fa";
			open SINGLETON, ">$folder/assem_id$assem_id.s.fa";

		}
	}

}

sub set_read {
	my ($l) = @_;

	my @a = split " ", $l;

	my $seq = "";
	if ($a[4] eq "D") {
		$seq = get_seq($a[2]);
	}
	else {
		$seq = rev_comp(get_seq($a[2]));
	}
	
	if ($a[9] == 1 or $a[9] == 3 or $a[9] == 4) {
		$first_reads{$a[3]} = $seq;
		$first_locations{$a[3]} = $a[1];
	}
	elsif ($a[9] == 2 or $a[9] == 6 or $a[9] == 7) {
		$second_reads{$a[3]} = $seq;
		$second_locations{$a[3]} = $a[1];
	}
	elsif ($a[9] == 0) {
		$singletons{$a[3]} = $seq;
		$singleton_locations{$a[3]} = $a[1];
	}
	elsif ($a[9] == 5) {
		$first_missing{$a[3]} = $seq;
		$first_locations{$a[3]} = $a[1];
	}
	elsif ($a[9] == 8) {
		$second_missing{$a[3]} = $seq;
		$second_locations{$a[3]} = $a[1];
	}

}

sub get_seq {
	my ($align) = @_;

	### Clean alignment from brackets and read-base
	my $seq = "";
	for (my $i = 0; $i<length($align); $i++) {
		if (substr($align, $i, 1) eq "[") {
			if (substr($align, $i+2, 1) ne "-") {
				$seq .= substr($align, $i+2, 1);
			}
			$i+=3;
		}
		else {
			$seq .= substr($align, $i, 1);
		}
	}

	return $seq;
}

sub rev_comp {
        my ($seq) = @_;
        my $new = reverse($seq);
        $new =~ tr/acgtACGT/tgcaTGCA/;

        return $new;
}


exit(0);


