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
#  Module: Analysis::Assembly::Validation::parseBlast.pl
#  Purpose:
#  In:
#  Out:
#


use Bio::SearchIO;

my $usage = "$0 blastoutputfile readfile targetfile\n";
my $file = shift or die $usage;
my $fasta = shift or die $usage;
my $ref = shift or die $usage;

my %SEQ = ();
my %SEQREF = ();

my $in = new Bio::SearchIO ( -format => 'blast', -file  => $file );

get_seq($fasta, \%SEQ);
get_seq($ref, \%SEQREF);

################################################################

open NOHIT, ">reads_nohit.fa";
open MISMATCH, ">reads_mismatch.fa";
open MISMATCH_SEVERE, ">reads_mismatch_severe.fa";
open UNRESOLVED, ">reads_unresolved.fa";
open MATCH, ">reads_match.fa";
open MATCH_GAPPED, ">reads_match_gapped.fa";
open CENT, ">reads_centormere.fa";
open ORG, ">reads_organelles.fa";

my $count = 0;
my $count_id = 0;
my $count_id_gapped = 0;
my $count_nohit = 0;
my $count_organelles = 0;
my $count_centromere = 0;
my $count_unresolved = 0;
my $count_mismatches = 0;
my $count_mismatches_severe = 0;

my $nuc_count_id = 0;
my $nuc_count_mis = 0;

RESULTS: while (my $blastResult = $in->next_result) {
	$count++;

	# What if no hit?
	if ($blastResult->hits == 0) {
		$count_nohit++;	
		print NOHIT ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
		next RESULTS;
	}

	######################################################
	# Analyze hit and set variables
	my $hit 	= $blastResult->next_hit;
	my $subjectName	= $hit->name;
	my $hsp 	= $hit->next_hsp;

	######################################################
	# Contamination hits:
	if ($subjectName eq "6" or $subjectName eq "7") {
		print ORG ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
		$count_organelles++;
		next RESULTS;
	}
	elsif (substr($subjectName, 0, 10) eq "Centromere") {
		print CENT ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
		$count_centromere++;
		next RESULTS;
	}

	######################################################
	# Correct hit:
	if ($hsp->num_identical == $hsp->hsp_length and $blastResult->query_length == $hsp->num_identical ) {
		print MATCH ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
		$count_id++;
		$nuc_count_id += $hsp->num_identical;
		next RESULTS;
	}

	######################################################
	# Get more information about the hit
	# count mismatches resulting from Ns in the sequences
	my $hsp_ambi_mismatches = get_ambi_mismatches($hsp->query_string, $hsp->hit_string, $hsp->homology_string);
	my $hsp_nuc_query = get_num_nuc($hsp->query_string);
	# ... and add the Ns at the beginning and end of the alignment that are not shown in the blast alignment
	my $hsp_ambi_up = 0;
	my $hsp_ambi_down = 0;
	if ($hsp->start('query') > 1) {
		my $r_length = min($hsp->start('hit')-1, $hsp->start('query')-1);
		my $seq_n = substr($SEQREF{$hit->name}, $hsp->start('hit')-$r_length-1, $r_length);
# print STDERR $seq_n, "\n";
		if (check_n($seq_n) == length($seq_n)) {
			$hsp_ambi_up += length($seq_n);
		}	
	}
	if ($hsp->end('query') < $blastResult->query_length) {
		my $r_length = min(length($SEQREF{$hit->name}) - $hsp->end('hit'), $blastResult->query_length - $hsp->end('query'));
		my $seq_n = substr($SEQREF{$hit->name}, $hsp->end('hit'), $r_length);
		if (check_n($seq_n) == length($seq_n)) {
			$hsp_ambi_down += length($seq_n);
		}
	}	

	######################################################
	# Correct hit with ambi mismatches:
	if (    $hsp->num_identical <= $hsp->hsp_length+$hsp_ambi_mismatches and
		$hsp->num_identical >= $hsp->hsp_length-$hsp_ambi_mismatches and
                min($hsp->num_identical+$hsp_ambi_mismatches, $hsp_nuc_query) + $hsp_ambi_up + $hsp_ambi_down == $blastResult->query_length) {
		print MATCH ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
		$count_id++;
		$nuc_count_id += $hsp->num_identical;
		next RESULTS;
	}

	######################################################
	# Correct hit but only nearly complete
# print STDERR $hsp->num_identical, "\t", $hsp_ambi_mismatches, "\t", $hsp_ambi_up, "\t", $hsp_ambi_down, "\t", $blastResult->query_length, "\n";
# print STDERR $hsp->length, "\n";
	if (	$hsp->num_identical+$hsp_ambi_mismatches+$hsp_ambi_up+$hsp_ambi_down >= $blastResult->query_length - 2 and
		$hsp->length <= $blastResult->query_length + 2) {
# print "FIRST HSP\n";
		print MISMATCH ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
		$count_mismatches++;
		$nuc_count_id += $hsp->num_identical;
		# Mismatches in the hsp + missing bases
 		$nuc_count_mis += ($hsp->hsp_length - $hsp->num_identical - $hsp_ambi_mismatches) + ($blastResult->query_length - $hsp_nuc_query);
		next RESULTS;
	}

	######################################################
	# more than minor difference but length somewhat represented
	# but length not aligned would not be enough for a sceond hit
	if (	$hsp->num_identical+$hsp_ambi_mismatches+$hsp_ambi_up+$hsp_ambi_down >= $blastResult->query_length - 20 and
		$hsp->length <= $blastResult->query_length + 20
	) {
		print MISMATCH_SEVERE ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
		$count_mismatches_severe++;
		next RESULTS;
	}

	######################################################
	# Catch the case where the first 2 hsp represent the hit but are splitted due to N stretches in one of the reads
	# Collect data
	my $hsp2 = "";
	my $hsp2_valid = 0;
	my $seq_n = "A"; # dummy
	my $ref_n = "A"; # dummy
	my $hsp2_ambi_mismatches = 0;
	my $hsp2_ambi_up = 0;
	my $hsp2_ambi_down = 0;

	my $hsp2_nuc_query = 0;

	######################################################
	# Does a meaningful hsp2 exist?
	if (	$hsp2 = $hit->next_hsp and
		$hsp->query->strand eq $hsp2->query->strand and # same orientation as hsp
		($hsp->end('query') < $hsp2->start('query') or $hsp->start('query') > $hsp2->end('query')) and 	# Neither read nor genome hit overlap but reads needs to be there completely
		($hsp->end('hit') < $hsp2->start('hit') or $hsp->start('hit') > $hsp2->end('hit')))
	{
		
		# get ambi count in hsp2
		$hsp2_ambi_mismatches = get_ambi_mismatches($hsp2->query_string, $hsp2->hit_string, $hsp2->homology_string);
		if ($hsp2->start('query') > 1) {
			my $tmp_seq_n = substr($SEQ{$blastResult->query_name}, 0, $hsp2->start('query')-1);
			if (check_n($tmp_seq_n) == length($tmp_seq_n)) {
				$hsp2_ambi_up += length($tmp_seq_n);
			}	
		}
		if ($hsp2->end('query') < $blastResult->query_length) {
			my $tmp_seq_n = substr($SEQ{$blastResult->query_name}, $hsp2->end('query'), $blastResult->query_length - $hsp2->start('query'));
			if (check_n($tmp_seq_n) == length($tmp_seq_n)) {
				$hsp2_ambi_down += length($tmp_seq_n);
			}
		}

		# Sequences in between are N in read or target
		$seq_n = get_n($SEQ{$blastResult->query_name}, $hsp->start('query'), $hsp->end('query'), $hsp2->start('query'), $hsp2->end('query'));
		$ref_n = get_n($SEQREF{$subjectName}, $hsp->start('hit'), $hsp->end('hit'), $hsp2->start('hit'), $hsp2->end('hit'));

        	$hsp2_nuc_query = get_num_nuc($hsp2->query_string);

# print $seq_n, "\n";
# print $ref_n, "\n";

		if (check_n($seq_n) == length($seq_n) or check_n($ref_n) == length($ref_n)) {
			$hsp2_valid = 1;
		}
	}

# print "HAS VALID:", $hsp2_valid,"\n";

	######################################################
	# Check combined alignments
	if ($hsp2_valid == 1) {
		my $all_ambi_up = max($hsp_ambi_up, $hsp2_ambi_up);
		my $all_ambi_down = max($hsp_ambi_down, $hsp2_ambi_down);

		# identical hit
		if ($hsp->num_identical+$hsp2->num_identical+($hsp_ambi_mismatches+$hsp2_ambi_mismatches+$all_ambi_up+$all_ambi_down)+length($seq_n) == $blastResult->query_length) {
			print MATCH_GAPPED ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
			$count_id_gapped++;
			$nuc_count_id += $hsp->num_identical+$hsp2->num_identical;
#print STDERR $seq_n, "\n";
#print STDERR $ref_n, "\n";
			next RESULTS;
		}

		# minor number of mismatches
		if (	$hsp->num_identical+$hsp2->num_identical+$hsp_ambi_mismatches+$hsp2_ambi_mismatches+$all_ambi_up+$all_ambi_down+length($seq_n) >= $blastResult->query_length-2) {
			print MISMATCH ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
			$count_mismatches++;
			$nuc_count_id += $hsp->num_identical+$hsp2->num_identical;
			$nuc_count_mis += $hsp->hsp_length - $hsp->num_identical - $hsp_ambi_mismatches + $hsp2->hsp_length - $hsp2->num_identical - $hsp2_ambi_mismatches;
                	# Mismatches in the hsp + missing bases
#print $blastResult->query_name, "\n";
#print $nuc_count_mis, "\t";
                	$nuc_count_mis += ($hsp->hsp_length - $hsp->num_identical - $hsp_ambi_mismatches);
                	$nuc_count_mis += ($hsp2->hsp_length - $hsp2->num_identical - $hsp2_ambi_mismatches);
                	$nuc_count_mis += ($blastResult->query_length - $hsp_nuc_query - $hsp2_nuc_query - $all_ambi_up - $all_ambi_down - length($seq_n));
#print $nuc_count_mis, "\n";
			next RESULTS;
		}

		# more than minor mismatches but length somewhat kept
		if (	$hsp->num_identical+$hsp2->num_identical+$hsp_ambi_mismatches+$hsp2_ambi_mismatches+$all_ambi_up+$all_ambi_down+length($seq_n) >= $blastResult->query_length-20) {
			print MISMATCH_SEVERE ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
			$count_mismatches_severe++;
			next RESULTS;
		}
	
	}


	######################################################
	# Read did not match so far, check only first hsp again
# 	if ($hsp->length('hit') >= $blastResult->query_length) {
# 		if ($hsp->num_identical >= $hsp->length('hit') -2) {
# 			print MISMATCH ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
# 			$count_mismatches++;
# 		}
# 		else {
# 			print MISMATCH_SEVERE ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
# 			$count_mismatches_severe++;
# 		}
# 	}
# 	else {
		print UNRESOLVED ">", $blastResult->query_name, "\n", $SEQ{$blastResult->query_name}, "\n";
		$count_unresolved++;
# 	}

}

print "=> count:\t\t\t", $count, "\n";
print "   count organellic:\t\t", $count_organelles, "\n";
print "   count centromere:\t\t", $count_centromere, "\n";
print "=> relevant count:\t\t", $count-($count_organelles+$count_centromere), "\t100%\n";
print "   count identical:\t\t", $count_id, "\t"; printf ("%.3f\n", 100*($count_id/($count-($count_organelles+$count_centromere))));
print "   count identical gapped:\t", $count_id_gapped, "\t"; printf ("%.3f\n", 100*($count_id_gapped/($count-($count_organelles+$count_centromere))));
print "   count mismatches:\t\t", $count_mismatches, "\t"; printf ("%.3f\n", 100*($count_mismatches/($count-($count_organelles+$count_centromere))));
print "   count mismatches(severe):\t", $count_mismatches_severe, "\t"; printf ("%.3f\n", 100*($count_mismatches_severe/($count-($count_organelles+$count_centromere))));
print "   count unresolved:\t\t", $count_unresolved, "\t"; printf ("%.3f\n", 100*($count_unresolved/($count-($count_organelles+$count_centromere))));
print "   count nohit:\t\t\t", $count_nohit, "\t"; printf ("%.3f\n", 100*($count_nohit/($count-($count_organelles+$count_centromere))));
print "########\n";
print "Identical, identical gapped and mismatching reads:\n";
print "   Nucleotides queried:\t", $nuc_count_id+$nuc_count_mis, "\t100%", "\n";
print "   Correct:\t\t", $nuc_count_id, "\t"; printf ("%.3f\n", 100*($nuc_count_id/($nuc_count_id+$nuc_count_mis))) ;
print "   Mismatch:\t\t", $nuc_count_mis, "\t"; printf ("%.3f\n", 100*($nuc_count_mis/($nuc_count_id+$nuc_count_mis)));


close NOHIT;
close MISMATCH;

sub get_ambi_mismatches {
	my ($query_string, $hit_string, $homology_string) = @_;

	my $ambiguous_mismatches = 0;
	for (my $p = 0; $p < length($query_string); $p++) {
		if (substr($homology_string, $p, 1) ne '\|') {
			if (substr($query_string, $p, 1) eq "N" or substr($hit_string, $p, 1) eq "N" or substr($query_string, $p, 1) eq "n" or substr($hit_string, $p, 1) eq "n") {
				$ambiguous_mismatches++;
			}
		}	
	}

	return $ambiguous_mismatches;
}

sub get_num_nuc {
	my ($s) = @_;
	$s =~ s/-//g;
	return length($s);
}

sub get_n {
	my ($s, $s1, $e1, $s2, $e2) = @_;
	my $seq = "";
	if ($e1 < $s2) {
		$seq = substr($s, $e1, $s2-$e1-1);
	}
	elsif ($e2 < $s1) {
		$seq = substr($s, $e2, $s1-$e2-1);
	}
	
	return $seq;
}

sub check_n {
	my ($seq) = @_;

	my $n = 0;
	for (my $i = 0; $i < length($seq); $i++) {
		if (substr($seq, $i, 1) eq "N" or substr($seq, $i, 1) eq "n") {
			$n++;
		}
	}
	
	return $n;
}

sub get_seq {
	my ($file, $hash_ref) = @_;
	my $id = "";
	my $seq = "";
	open FASTA, $file or die $usage;
	while (my $line = <FASTA>) {
        	chomp($line);
	        if (substr($line, 0, 1) eq ">") {
        	        if ($seq ne "") {
                	        ${$hash_ref}{$id} = $seq;
	                }
        	        $seq = "";
                	my @a = split " ", $line;
	                $id = substr($a[0], 1, length($a[0])-1);
	        }
        	else {
                	$seq .= $line;
	        }
	}
	close FASTA;
	${$hash_ref}{$id} = $seq;
}

sub max {
	my ($a, $b) = @_;

	return $a if $a > $b;
	return $b;
}

sub min {
	my ($a, $b) = @_;

	return $a if $a < $b;
	return $b;
}

