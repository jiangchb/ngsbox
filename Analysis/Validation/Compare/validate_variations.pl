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
#  Module: Analysis::Validation::Compare::validate_variations.pl
#  Purpose:
#  In:
#  Out:
#


use lib "$ENV{PGSP}/Prediction/Validation/";

use Genome;
use Compare;
use VariationListing;


open SNP1, ">overlapping.snp1.txt";
open SNP2, ">overlapping.snp2.txt";

my $usage = "$0 refseq snpfile1 deletionfile1 insertionfile1 snpfile2 deletionfile2 insertionfile2\n";

my $reffile = shift or die $usage;
my $snp1file = shift or die $usage;
my $deletion1file = shift or die $usage;
my $insertion1file = shift or die $usage;
my $snp2file = shift or die $usage;
my $deletion2file = shift or die $usage;
my $insertion2file = shift or die $usage;

my $genome1 = new Genome();
my $genome2 = new Genome();
print STDERR "Get genome1...\n";
$genome1->get_listings($reffile, $snp1file, $deletion1file, $insertion1file);
print STDERR "...done\n";
print STDERR "Get genome2...\n";
$genome2->get_listings($reffile, $snp2file, $deletion2file, $insertion2file);
print STDERR "...done\n";

my $variation1 = new VariationListing();
my $variation2 = new VariationListing();
print STDERR "Get variation list 1...\n";
$variation1->calc_private_variation_sets($snp1file, $deletion1file, $insertion1file, $snp2file, $deletion2file, $insertion2file, 0);
print STDERR "...done\n";
print STDERR "Get variation list 2...\n";
$variation2->calc_private_variation_sets($snp2file, $deletion2file, $insertion2file, $snp1file, $deletion1file, $insertion1file, 1);
print STDERR "...done\n";

my $compare = new Compare();

# Parse variation list
my $count_snp1_priv_ident = 0;
my $count_snp1_priv_diff = 0;
my $count_del1_priv_ident = 0;
my $count_del1_priv_diff = 0;
my $count_ins1_priv_ident = 0;
my $count_ins1_priv_diff = 0;
my $count_del1_priv_ident_bp = 0;
my $count_del1_priv_diff_bp = 0;
my $count_ins1_priv_ident_bp = 0;
my $count_ins1_priv_diff_bp = 0;

my $count_del1_priv_ident_1bp = 0;
my $count_del1_priv_diff_1bp = 0;
my $count_ins1_priv_ident_1bp = 0;
my $count_ins1_priv_diff_1bp = 0;
my $count_del1_priv_ident_4bp = 0;
my $count_del1_priv_diff_4bp = 0;
my $count_ins1_priv_ident_4bp = 0;
my $count_ins1_priv_diff_4bp = 0;

my $count_snp2_priv_ident = 0;
my $count_snp2_priv_diff = 0;
my $count_del2_priv_ident = 0;
my $count_del2_priv_diff = 0;
my $count_ins2_priv_ident = 0;
my $count_ins2_priv_diff = 0;
my $count_del2_priv_ident_bp = 0;
my $count_del2_priv_diff_bp = 0;
my $count_ins2_priv_ident_bp = 0;
my $count_ins2_priv_diff_bp = 0;

my $count_del2_priv_ident_1bp = 0;
my $count_del2_priv_diff_1bp = 0;
my $count_ins2_priv_ident_1bp = 0;
my $count_ins2_priv_diff_1bp = 0;
my $count_del2_priv_ident_4bp = 0;
my $count_del2_priv_diff_4bp = 0;
my $count_ins2_priv_ident_4bp = 0;
my $count_ins2_priv_diff_4bp = 0;


my $min_offset = 5;
my $max_offset = 36;

foreach my $key (keys ( %{$variation1->{snp}})) {
	my ($chr, $pos) = split "#", $key;
	my $offset = $min_offset;
	my $comp = 0;
	while ($comp == 0 and $offset < $max_offset) {
		$offset++;
		my $seq1 = $genome1->get_enlarged_sequence_from_listing($chr, $pos, $pos, $offset, $pos, $pos);
		my $seq2 = $genome2->get_enlarged_sequence_from_listing($chr, $pos, $pos, $offset, $pos, $pos);
		$comp = $genome2->compare($seq1, $seq2, $chr, $pos-($offset), $pos+($offset));
	}
	if ($comp == 1) { $count_snp1_priv_ident++; print SNP1 $chr, "\t", $pos, "\n"; } 
	else 		{ $count_snp1_priv_diff++; }
}

foreach my $key (keys ( %{$variation1->{del}})) {
        my ($chr, $begin, $end) = split "#", $key;
	my $offset = $min_offset;
        my $comp = 0;
        while ($comp == 0 and $offset < $max_offset) {
		$offset++;
	        my $seq1 = $genome1->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
        	my $seq2 = $genome2->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
	        $comp = $genome2->compare($seq1, $seq2, $chr, $begin-($offset), $end+($offset));
	}
	my $len = $variation1->{del}{$key};
        if ($comp == 1) { 
		$count_del1_priv_ident++; 
		$count_del1_priv_ident_bp+=$len; 
		if ($len <= 3) 	{ $count_del1_priv_ident_1bp++; }
		else 		{ $count_del1_priv_ident_4bp++; }
	}
        else { 
		$count_del1_priv_diff++; 
		$count_del1_priv_diff_bp+=$len;
		if ($len <= 3) 	{ $count_del1_priv_diff_1bp++; }
                else		{ $count_del1_priv_diff_4bp++; }
	}
}

foreach my $key (keys ( %{$variation1->{ins}})) {
        my ($chr, $begin, $end) = split "#", $key;
	my $offset = $min_offset;
        my $comp = 0;
        while ($comp == 0 and $offset < $max_offset) {
                $offset++;
	        my $seq1 = $genome1->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
        	my $seq2 = $genome2->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
	        $comp = $genome2->compare($seq1, $seq2, $chr, $begin-($offset), $end+($offset));
	}
	my $len = length($variation1->{ins}{$key});
        if ($comp == 1) { 
		$count_ins1_priv_ident++; 
		$count_ins1_priv_ident_bp+=$len;
		if ($len <= 3) 	{ $count_ins1_priv_ident_1bp++; }
                else 		{ $count_ins1_priv_ident_4bp++; }
	}
        else { 
		$count_ins1_priv_diff++; 
		$count_ins1_priv_diff_bp+=$len;
		if ($len <= 3) 	{ $count_ins1_priv_diff_1bp++; }
                else 		{ $count_ins1_priv_diff_4bp++; }
	}
}

foreach my $key (keys ( %{$variation2->{snp}})) {
        my ($chr, $pos) = split "#", $key;
	my $offset = $min_offset;
        my $comp = 0;
        while ($comp == 0 and $offset < $max_offset) {
                $offset++;
	        my $seq1 = $genome2->get_enlarged_sequence_from_listing($chr, $pos, $pos, $offset, $pos, $pos);
        	my $seq2 = $genome1->get_enlarged_sequence_from_listing($chr, $pos, $pos, $offset, $pos, $pos);
	        $comp = $genome1->compare($seq1, $seq2, $chr, $pos-($offset), $pos+($offset));
	}
        if ($comp == 1) { $count_snp2_priv_ident++; print SNP2 $chr, "\t", $pos, "\n"; }
        else            { $count_snp2_priv_diff++; }
}

foreach my $key (keys ( %{$variation2->{del}})) {
        my ($chr, $begin, $end) = split "#", $key;
	my $offset = $min_offset;
        my $comp = 0;
        while ($comp == 0 and $offset < $max_offset) {
                $offset++;
	        my $seq1 = $genome2->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
        	my $seq2 = $genome1->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
	        $comp = $genome1->compare($seq1, $seq2, $chr, $begin-($offset), $end+($offset));
	}
	my $len = $variation2->{del}{$key};
        if ($comp == 1) { 
		$count_del2_priv_ident++; 
		$count_del2_priv_ident_bp+=$len; 
		if ($len <=3) 	{ $count_del2_priv_ident_1bp++; }
                else 		{ $count_del2_priv_ident_4bp++; }
	}
        else { 
		$count_del2_priv_diff++; 
		$count_del2_priv_diff_bp+=$len;
		if ($len <= 3) 	{ $count_del2_priv_diff_1bp++; }
                else		{ $count_del2_priv_diff_4bp++; } 
	} 
}

foreach my $key (keys ( %{$variation2->{ins}})) {
        my ($chr, $begin, $end) = split "#", $key;
	my $offset = $min_offset;
        my $comp = 0;
        while ($comp == 0 and $offset < $max_offset) {
                $offset++;
        	my $seq1 = $genome2->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
	        my $seq2 = $genome1->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
        	$comp = $genome1->compare($seq1, $seq2, $chr, $begin-($offset), $end+($offset));
	}
	my $len = length($variation2->{ins}{$key});
        if ($comp == 1) { 
		$count_ins2_priv_ident++; 
		$count_ins2_priv_ident_bp+=$len; 
		if ($len <= 3) 	{ $count_ins2_priv_ident_1bp++; }
                else 		{ $count_ins2_priv_ident_4bp++; }
	}
        else { 
		$count_ins2_priv_diff++; 
		$count_ins2_priv_diff_bp+=$len;
		if ($len <= 3) 	{ $count_ins2_priv_diff_1bp++; }
                else 		{ $count_ins2_priv_diff_4bp++; }
	}
}


# Print
print "Set1:\n";
print "  SNPs:\n";
print "    Total:\t\t", $variation1->{snp_num}, "\n";
print "    Identical to Set1:\t", $variation1->{snp_num} - $variation1->{snp_num_priv}, "\n";
print "    Same in context:\t", $count_snp1_priv_ident, "\n";
print "    Different:\t\t", $count_snp1_priv_diff, "\n\n";

print "  Dels:\n";
print "    Total:\t\t", $variation1->{del_num}, "\n";
print "    Identical to Set1:\t", $variation1->{del_num} - $variation1->{del_num_priv}, "\n";
print "    Same in context:\t", $count_del1_priv_ident, "\t", $count_del1_priv_ident_bp, "\n";
print "    Different:\t\t", $count_del1_priv_diff, "\t", $count_del1_priv_diff_bp, "\n\n";

print "    Identical 1-3bp:\t", $variation1->{del_num_id_1bp}, "\n";
print "    Same in context 1-3bp:\t", $count_del1_priv_ident_1bp, "\n";
print "    Different 1-3bp:\t\t", $count_del1_priv_diff_1bp, "\n\n";

print "    Identical >= 4bp:\t", $variation1->{del_num_id_4bp}, "\n";
print "    Same in context >= 4bp:\t", $count_del1_priv_ident_4bp, "\n";
print "    Different >= 4bp:\t\t", $count_del1_priv_diff_4bp, "\n\n";

print "  Ins:\n";
print "    Total:\t\t", $variation1->{ins_num}, "\n";
print "    Identical to Set1:\t", $variation1->{ins_num} - $variation1->{ins_num_priv}, "\n";
print "    Same in context:\t", $count_ins1_priv_ident, "\t", $count_ins1_priv_ident_bp, "\n";
print "    Different:\t\t", $count_ins1_priv_diff, "\t", $count_ins1_priv_diff_bp, "\n\n";

print "    Identical 1-3bp:\t", $variation1->{ins_num_id_1bp}, "\n";
print "    Same in context 1-3bp:\t", $count_ins1_priv_ident_1bp, "\n";
print "    Different 1-3bp:\t\t", $count_ins1_priv_diff_1bp, "\n\n";

print "    Identical >= 4bp:\t", $variation1->{ins_num_id_4bp}, "\n";
print "    Same in context >= 4bp:\t", $count_ins1_priv_ident_4bp, "\n";
print "    Different >= 4bp:\t\t", $count_ins1_priv_diff_4bp, "\n\n";



print "Set2:\n";
print "  SNPs:\n";
print "    Total:\t\t", $variation2->{snp_num}, "\n";
print "    Identical to Set1:\t", $variation2->{snp_num} - $variation2->{snp_num_priv}, "\n";
print "    Same in context:\t", $count_snp2_priv_ident, "\n";
print "    Different:\t\t", $count_snp2_priv_diff, "\n\n";

print "  Dels:\n";
print "    Total:\t\t", $variation2->{del_num}, "\n";
print "    Identical to Set1:\t", $variation2->{del_num} - $variation2->{del_num_priv}, "\n";
print "    Same in context:\t", $count_del2_priv_ident, "\t", $count_del2_priv_ident_bp, "\n";
print "    Different:\t\t", $count_del2_priv_diff, "\t", $count_del2_priv_diff_bp, "\n\n";

print "    Identical 1-3bp:\t", $variation2->{del_num_id_1bp}, "\n";
print "    Same in context 1-3bp:\t", $count_del2_priv_ident_1bp, "\n";
print "    Different 1-3bp:\t\t", $count_del2_priv_diff_1bp, "\n\n";

print "    Identical >= 4bp:\t", $variation2->{del_num_id_4bp}, "\n";
print "    Same in context >= 4bp:\t", $count_del2_priv_ident_4bp, "\n";
print "    Different >= 4bp:\t\t", $count_del2_priv_diff_4bp, "\n\n";

print "  Ins:\n";
print "    Total:\t\t", $variation2->{ins_num}, "\n";
print "    Identical to Set1:\t", $variation2->{ins_num} - $variation2->{ins_num_priv}, "\n";
print "    Same in context:\t", $count_ins2_priv_ident, "\t", $count_ins2_priv_ident_bp, "\n";
print "    Different:\t\t", $count_ins2_priv_diff, "\t", $count_ins2_priv_diff_bp, "\n\n";

print "    Identical 1-3bp:\t", $variation2->{ins_num_id_1bp}, "\n";
print "    Same in context 1-3bp:\t", $count_ins2_priv_ident_1bp, "\n";
print "    Different 1-3bp:\t\t", $count_ins2_priv_diff_1bp, "\n\n";

print "    Identical >= 4bp:\t", $variation2->{ins_num_id_4bp}, "\n";
print "    Same in context >= 4bp:\t", $count_ins2_priv_ident_4bp, "\n";
print "    Different >= 4bp:\t\t", $count_ins2_priv_diff_4bp, "\n\n";



