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
#  Module: Analysis::Validation::Compare::Genome.pm
#  Purpose:
#  In:
#  Out:
#



package Genome;

sub new {
	my $self = {
		genome => {},	
		chr_size => {},
		chr_seq => {},
		snp_list => {},
		del_list => {},
		ins_list => {},
	};
	bless $self;
	return $self;
}

#############################################################################################
### Get chr seq
sub get_listings {
	my ($self, $reffile, $snpfile, $delfile, $insfile) = @_;

	open FILE, $reffile or die "Cannot open file $reffile\n";
	my $curr_chr = "";
	while (my $line = <FILE>) {
		chomp($line);
		if (substr($line, 0, 1) eq ">") {
                        if ($curr_chr ne "") {
				print STDERR "Got Chromosome ", $curr_chr, " size: ", $self->{chr_size}{$curr_chr}, "\n";
			}
			my @a = split " ", $line;
                        $curr_chr = substr($a[0], 1, length($a[0])-1);
			$self->{chr_size}{$curr_chr} = 0;
                }
		else {
                        if ($curr_chr ne "") {
				$self->{chr_seq}{$curr_chr} .= $line;
				$self->{chr_size}{$curr_chr} += length($line);
                        }
                }
	}
	if ($curr_chr ne "") {
		print STDERR "Got Chromosome ", $curr_chr, " size: ", $self->{chr_size}{$curr_chr}, "\n";
        }

	# SNP file
	open FILE, $snpfile or die "Cannot open file $snpfile\n";
	while (my $line = <FILE>) {
                my @a = split " ", $line;
		$self->{snp_list}{$a[1]}{$a[2]} = $a[4];
	}
	close FILE;

	# Deletion file
	open FILE, $delfile or die "Cannot open file $delfile\n";
	while (my $line = <FILE>) {
                my @a = split " ", $line;
		for (my $i = $a[2]; $i <= $a[3]; $i++) {
	                $self->{del_list}{$a[1]}{$i} = $a[1]."#".$a[2]."#".$a[3];
		}
        }
        close FILE;

	# Insertion file
	open FILE, $insfile or die "Cannot open file $insfile\n";
        while (my $line = <FILE>) {
                my @a = split " ", $line;
		$self->{ins_list}{$a[1]}{$a[2]} .= $a[5];
        }
        close FILE;


}


#########################################################################################
###
sub get_enlarged_sequence_from_listing {
	my ($self, $chr, $begin, $end_var, $offset, $begin_check, $end_var_check) = @_;
	
	my $seq = "";

	# Set start position
	my $start = $begin;
 	my $homopolymer = 1;
 	START: while ($homopolymer == 1) {
 		if ($start < 1) {
 			$homopolymer = 0;
 		}
 		else {
			if (substr($self->{chr_seq}{$chr}, $start-2, 1) ne substr($self->{chr_seq}{$chr}, $start-1, 1)) {
				$homopolymer = 0;
			}
 		}
 		$start--;
 	}
	if ($begin-$offset < $start) {
		$start = $begin-$offset;
	}
	if ($start < 1) {
		$start = 1;
	}

	# Set end position
	my $end = $end_var;
	$homopolymer = 1;
	ENDE: while ($homopolymer == 1) {
		if (not defined($self->{chr_size}{$chr})) {
			print STDERR "Cannot find the chromsome idenfier: \"", $chr, "\"\n";
			print STDERR "Known identifiers are:\n";
			foreach my $key (keys %{$self->{chr_size}}) {
				print STDERR $key, "\n";
			}
			exit(1);
		}
                if ($end > $self->{chr_size}{$chr}) {
                         $homopolymer = 0;
                }
                else {
			if (substr($self->{chr_seq}{$chr}, $end-1, 1) ne substr($self->{chr_seq}{$chr}, $end, 1)) {
				$homopolymer = 0;
			}
                }
                $end++;
	}
	if ($end_var+$offset > $end) {
		$end = $end_var + $offset;
	}
	if ($end > $self->{chr_size}{$chr} - 1) {
		$end = $self->{chr_size}{$chr} - 1;
	}
	

	# Get sequence
	for (my $i = $start; $i <= $end; $i++) {
		# Is this position deleted?
		if (defined($self->{del_list}{$chr}{$i})) {
		}
		else {
			if (defined($self->{snp_list}{$chr}{$i})) {
				if ($i < $begin_check or $i > $end_var_check) {
					$seq .= "N";
				}
				else {
					$seq .= $self->{snp_list}{$chr}{$i};
				}
			}
			else {
				$seq .= substr($self->{chr_seq}{$chr}, $i-1, 1);
			}

			# Is there a insertion				
			if (defined($self->{ins_list}{$chr}{$i})) {
				$seq .= $self->{ins_list}{$chr}{$i};
			}
			
		}
	}

	return $seq;
}


sub compare {
	my ($self, $seq1, $seq2, $chr, $begin, $end) = @_;
	
	# Check if there is a change at all
	my $no_change_flag = 1;
	for (my $i = $begin; $i <= $end; $i++) {
		if (defined($self->{ins_list}{$chr}{$i})) {
			$no_change_flag = 0;
		}
		if (defined($self->{snp_list}{$chr}{$i})) {
                        $no_change_flag = 0;
                }
		if (defined($self->{del_list}{$chr}{$i})) {
                        $no_change_flag = 0;
                }
	}

	if ($no_change_flag == 1) {
		return 0;
	}

	# Check seqeunces
        my $len_diff = abs(length($seq1) - length($seq2));

        for (my $i = 0; $i <= $len_diff; $i++) {
                my $start1 = 0;
                my $start2 = 0;
                my $len;
                if (length($seq1) > length($seq2)) {
                        $start1 = $i;
                        $len = length($seq2);
                }
                elsif (length($seq1) < length($seq2)) {
                        $start2 = $i;
                        $len = length($seq1);
                }
                else {
                        $len = length($seq1);
                }

                my $ident = 1;

                for (my $j = 0; $j < $len; $j++) {
                        if (    substr($seq1, $start1+$j, 1) ne substr($seq2, $start2+$j, 1) and
                                substr($seq1, $start1+$j, 1) ne "N" and
                                substr($seq2, $start2+$j, 1) ne "N"
                        ) {
                                $ident = 0;
                                $j = $len;
                        }
                }

                if ($ident == 1) {
                        return 1;
                }
        }

        return 0;
}



#############################################################################################
#############################################################################################
#############################################################################################
## This is based on a pre-calculation of the genomes
## Which does not allow for the flexible adjustment
## of variations nearby.
#############################################################################################
### Calculates the genome sequence of a sampled genome
sub calc_chr_seq {
	my ($self, $reffile, $snpfile, $delfile, $insfile) = @_;

	# Reference sequence
	open FILE, $reffile or die "Cannot open file $reffile\n";
	my $curr_chr = "";
	my $curr_pos = 0;
	while (my $line = <FILE>) {
		chomp($line);
		if (substr($line, 0, 1) eq ">") {
			if ($curr_chr ne "") {
				$self->{chr_size}{$curr_chr} = $curr_pos;
				print STDERR "Got Chromosome ", $curr_chr, " size: ", $self->{chr_size}{$curr_chr}, "\n";;
			}
			my @a = split " ", $line;
			$curr_chr = substr($a[0], 1, length($a[0])-1);
			$curr_pos = 0;
		}
		else {
			if ($curr_chr ne "") {
				for (my $i = 0; $i < length($line); $i++) {
					$curr_pos++;
					$self->{genome}{$curr_chr}[$curr_pos] = substr($line, $i, 1);
				}
			}		
		}
	}
	if ($curr_chr ne "") {
        	$self->{chr_size}{$curr_chr} = $curr_pos;
		print STDERR "Got Chromosome ", $curr_chr, " size: ", $self->{chr_size}{$curr_chr}, "\n";
        }
	close FILE;

	# SNP file
	open FILE, $snpfile or die "Cannot open file $snpfile\n";
	while (my $line = <FILE>) {
                my @a = split " ", $line;
		$self->{genome}{$a[1]}[$a[2]] = $a[4];
	}
	close FILE;

	# Deletion file
	open FILE, $delfile or die "Cannot open file $delfile\n";
	while (my $line = <FILE>) {
                my @a = split " ", $line;
		for (my $i = $a[2]; $i <= $a[3]; $i++) {
	                $self->{genome}{$a[1]}[$i] = "-";
		}
        }
        close FILE;

	# Insertion file
	open FILE, $insfile or die "Cannot open file $insfile\n";
        while (my $line = <FILE>) {
                my @a = split " ", $line;
		$self->{genome}{$a[1]}[$a[2]] .= $a[5];
        }
        close FILE;


}

sub get_sequence {
	my ($self, $chr, $begin, $end) = @_;

	my $seq = "";

	if (not defined($self->{chr_size}{$chr})) {
		print STDERR "Cannot find the chromsome idenfier: \"", $chr, "\"\n";
                        print STDERR "Known identifiers are:\n";
                     	foreach my $key (keys %{$self->{chr_size}}) {
                                print STDERR $key, "\n";
                        }
                        exit(1);
        }
        if ($end > $self->{chr_size}{$chr}) {
        	$end = $self->{chr_size}{$chr};
	}
	

	for (my $i = $begin; $i <= $end; $i++) {
                if ($self->{genome}{$chr}[$i] ne "-") {
                        $seq .= $self->{genome}{$chr}[$i];
                }
        }

	return $seq;
}
	
#########################################################################################
### Calculates CDS and protein changes for all SNPs in "SNP_lists" and fills SNP object
sub get_enlarged_sequence {
	my ($self, $chr, $begin, $end_var, $offset) = @_;
	
	my $seq = "";

	# Set start position
	my $start = $begin;
 	my $homopolymer = 1;
 	START: while ($homopolymer == 1) {
 		if ($start < 1) {
 			$homopolymer = 0;
 		}
 		else {
			if (substr($self->{genome}{$chr}[$start-1], 0, 1) ne substr($self->{genome}{$chr}[$start], 0, 1)) {
				$homopolymer = 0;
			}
 		}
 		$start--;
 	}
	if ($begin-$offset-1 < $start) {
		$start = $begin-$offset-1;
	}
	if ($start < 1) {
		$start = 1;
	}

	# Set end position
	my $end = $end_var;
	$homopolymer = 1;
	ENDE: while ($homopolymer == 1) {
		if (not defined($self->{chr_size}{$chr})) {
			print STDERR "Cannot find the chromsome idenfier: \"", $chr, "\"\n";
			print STDERR "Known identifiers are:\n";
			foreach my $key (keys %{$self->{chr_size}}) {
				print STDERR $key, "\n";
			}
			exit(1);
		}
                if ($end > $self->{chr_size}{$chr}) {
                         $homopolymer = 0;
                }
                else {
			if (substr($self->{genome}{$chr}[$end+1], 0, 1) ne substr($self->{genome}{$chr}[$end], 0, 1)) {
				$homopolymer = 0;
			}
                }
                $end++;
	}
	if ($end_var+$offset+1 > $end) {
		$end = $end_var + $offset + 1;
	}
	if ($end > $self->{chr_size}{$chr}) {
		$end = $self->{chr_size}{$chr};
	}
	

	# Get sequence
	for (my $i = $start; $i <= $end; $i++) {
		if ($self->{genome}{$chr}[$i] ne "-") {
			$seq .= $self->{genome}{$chr}[$i];
			#print $i, "\t", $seq, "\n";
		}
		else {
			#print "deletion\n";
		}
	}

	return $seq;
}

1;



