#! /usr/bin/perl
use strict;

use lib "$ENV{PGSP}/Prediction/Validation/";

my $min_offset = 5;
my $max_offset = 36;

use Genome;
use Compare;
use VariationListing;

my $usage = "Usage:\n$0\nrefseq.txt mn_fragments.txt mn_snps.txt mn_deletions.txt mn_insertions.txt snps.txt deletion.txt insertion.txt\n\n";
my $ref = shift or die $usage;
my $frag = shift or die $usage;
my $mnsnps = shift or die $usage;
my $mndel = shift or die $usage;
my $mnins = shift or die $usage;
my $snps = shift or die $usage;
my $del = shift or die $usage;
my $ins = shift or die $usage;

## Read in fragments		
open FILE, $frag or die $usage;
my %MNPOS = ();
my $mn_frag_trimming = 0;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	for (my $i = $a[1]+$mn_frag_trimming; $i <= $a[2]-$mn_frag_trimming; $i++) {
		$MNPOS{$a[0]."#".$i} = 1;
	}
}
close FILE;

# Set up genomes
print STDERR "Get listings...\n";
my $mngenome = new Genome();
#$mngenome->calc_chr_seq($ref, $mnsnps, $mndel, $mnins);
$mngenome->get_listings($ref, $mnsnps, $mndel, $mnins);
print STDERR "..done...\n";
my $genome = new Genome();
#$genome->calc_chr_seq($ref, $snps, $del, $ins);
$genome->get_listings($ref, $snps, $del, $ins);
print STDERR "..done\n";

my $mnvariation = new VariationListing();
$mnvariation->calc_private_variation_sets($mnsnps, $mndel, $mnins, $snps, $del, $ins);
my $variation = new VariationListing();
$variation->calc_private_variation_sets($snps, $del, $ins, $mnsnps, $mndel, $mnins);

my $compare = new Compare();

open SNPFP, ">snps.fp.txt";
open DELFP, ">del.fp.txt";
open INSFP, ">ins.fp.txt";

open SNPTP, ">snps.tp.txt";
open DELTP, ">del.tp.txt";
open INSTP, ">ins.tp.txt";

open SNPFN, ">snps.fn.txt";
open DELFN, ">del.fn.txt";
open INSFN, ">ins.fn.txt";


my $snps_tp = 0;
my $snps_fp = 0;
my $snps_fn = 0;

my $del_tp = 0;
my $del_fp = 0;
my $del_fn = 0;
my $del_tp_4bp = 0;
my $del_fp_4bp = 0;
my $del_fn_4bp = 0;


my $ins_tp = 0;
my $ins_fp = 0;
my $ins_fn = 0;
my $ins_tp_4bp = 0;
my $ins_fp_4bp = 0;
my $ins_fn_4bp = 0;


# Compare SNPs to MN
foreach my $key (keys ( %{$variation->{snp}})) {
	my ($chr, $pos) = split "#", $key;
	my $offset = $min_offset;
	my $comp = 0;
	if (defined($MNPOS{$chr."#".$pos})) {
	        while ($comp == 0 and $offset < $max_offset) {
        	        $offset++;
                	#my $seq1 = $genome->get_enlarged_sequence($chr, $pos, $pos, $offset);
	                #my $seq2 = $mngenome->get_enlarged_sequence($chr, $pos, $pos, $offset);
	               	my $seq1 = $genome->get_enlarged_sequence_from_listing($chr, $pos, $pos, $offset, $pos, $pos);
			my $seq2 = $mngenome->get_enlarged_sequence_from_listing($chr, $pos, $pos, $offset, $pos, $pos);
        	        #$comp = $compare->compare($seq1, $seq2);
			$comp = $genome->compare($seq1, $seq2, $chr, $pos-($offset), $pos+($offset));
	        }
		if ($comp != 1) { 
			$snps_fp++; 
			print SNPFP $chr, "\t", $pos, "\tFP\n";  
		}
	}
}

# Compare MN to SNPs
foreach my $key (keys ( %{$mnvariation->{snp}})) {
        my ($chr, $pos) = split "#", $key;
        my $offset = $min_offset;
        my $comp = 0;
        if (defined($MNPOS{$chr."#".$pos})) {
                while ($comp == 0 and $offset < $max_offset) {
                        $offset++;
                        #my $seq1 = $mngenome->get_enlarged_sequence($chr, $pos, $pos, $offset);
                        #my $seq2 = $genome->get_enlarged_sequence($chr, $pos, $pos, $offset);
			my $seq1 = $mngenome->get_enlarged_sequence_from_listing($chr, $pos, $pos, $offset, $pos, $pos);
                        my $seq2 = $genome->get_enlarged_sequence_from_listing($chr, $pos, $pos, $offset, $pos, $pos);
                        #$comp = $compare->compare($seq1, $seq2);
                        $comp = $genome->compare($seq1, $seq2, $chr, $pos-($offset), $pos+($offset));
                }
		if ($comp == 1) {
			$snps_tp++;
			print SNPTP $chr, "\t", $pos, "\tTP\n";
		}
                else  { 
			print SNPFN $chr, "\t", $pos, "\tFN\n";
			$snps_fn++; 
		}
        }
}

# Compare deletions
foreach my $key (keys ( %{$variation->{del}})) {
        my ($chr, $begin, $end) = split "#", $key;
        my $offset = $min_offset;
        my $comp = 0;
	if (defined($MNPOS{$chr."#".$begin}) and defined($MNPOS{$chr."#".$end})) {
		my $len = $variation->{del}{$key};
	        while ($comp == 0 and $offset < $max_offset) {
        	        $offset++;
                	#my $seq1 = $genome->get_enlarged_sequence($chr, $begin, $end, $offset);
	                #my $seq2 = $mngenome->get_enlarged_sequence($chr, $begin, $end, $offset);
			my $seq1 = $genome->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
                        my $seq2 = $mngenome->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
	       	        #$comp = $compare->compare($seq1, $seq2), "\n";
	       	        $comp = $genome->compare($seq1, $seq2, $chr, $begin-($offset), $end+($offset));
	        }
	        if ($comp != 1) { 
			$del_fp++; 
			print DELFP $chr, "\t", $begin, "\t", $end, "\tFP\t", $len, "\n"; 
			if ($len >= 4) {
				$del_fp_4bp++;
			}
		}
	}
}

foreach my $key (keys ( %{$mnvariation->{del}})) {
        my ($chr, $begin, $end) = split "#", $key;
        my $offset = $min_offset;
        my $comp = 0;
        if (defined($MNPOS{$chr."#".$begin}) and defined($MNPOS{$chr."#".$end})) {
		my $seq1 = ""; my $seq2 = "";
		my $len = $mnvariation->{del}{$key};
                while ($comp == 0 and $offset < $max_offset) {
                        $offset++;
                        #$seq1 = $genome->get_enlarged_sequence($chr, $begin, $end, $offset);
                        #$seq2 = $mngenome->get_enlarged_sequence($chr, $begin, $end, $offset);
			$seq1 = $genome->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
                        $seq2 = $mngenome->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
                        #$comp = $compare->compare($seq1, $seq2), "\n";
                        $comp = $genome->compare($seq1, $seq2, $chr, $begin-($offset), $end+($offset));
                }
		if ($comp == 1) {
                        $del_tp++;
                        if ($len >= 4) {
                                $del_tp_4bp++;
                        }
			print DELTP $chr, "\t", $begin, "\t", $end, "\tTP\t", $len, "\n";
                }
                else { 
			$del_fn++; 
			if ($len >= 4) {
				$del_fn_4bp++;
			}
			print DELFN $chr, "\t", $begin, "\t", $end, "\tFN\t", $len, "\n";
		}
        }
}

# Compare insertions
foreach my $key (keys ( %{$variation->{ins}})) {
        my ($chr, $begin, $end) = split "#", $key;
        my $offset = $min_offset;
        my $comp = 0;
        if (defined($MNPOS{$chr."#".$begin}) and defined($MNPOS{$chr."#".$end})) {
		my $len = length($variation->{ins}{$key});
                while ($comp == 0 and $offset < $max_offset) {
                        $offset++;
                        #my $seq1 = $genome->get_enlarged_sequence($chr, $begin, $end, $offset);
                        my $seq1 = $genome->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
                        #my $seq2 = $mngenome->get_enlarged_sequence($chr, $begin, $end, $offset);
                        my $seq2 = $mngenome->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
	                #$comp = $compare->compare($seq1, $seq2), "\n";
	                $comp = $genome->compare($seq1, $seq2, $chr, $begin-($offset), $end+($offset));
                }
                if ($comp != 1) { 
			$ins_fp++; 
			print INSFP $chr, "\t", $begin, "\t", $end, "\tFP\t", $len, "\n";
			if ($len >= 4) {
				$ins_fp_4bp++;
			}
		}
        }
}

foreach my $key (keys ( %{$mnvariation->{ins}})) {
        my ($chr, $begin, $end) = split "#", $key;
        my $offset = $min_offset;
        my $comp = 0;
        if (defined($MNPOS{$chr."#".$begin}) and defined($MNPOS{$chr."#".$end})) {
		my $seq1 = "";
		my $seq2 = "";
		my $len = length($mnvariation->{ins}{$key});
                while ($comp == 0 and $offset < $max_offset) {
                        $offset++;
                        #$seq1 = $genome->get_enlarged_sequence($chr, $begin, $end, $offset);
                        $seq1 = $genome->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
                        #$seq2 = $mngenome->get_enlarged_sequence($chr, $begin, $end, $offset);
                        $seq2 = $mngenome->get_enlarged_sequence_from_listing($chr, $begin, $end, $offset, $begin, $end);
                        #$comp = $compare->compare($seq1, $seq2), "\n";
                        $comp = $genome->compare($seq1, $seq2, $chr, $begin-($offset), $end+($offset));
                }
		if ($comp == 1) {
                        $ins_tp++;
                        if ($len >= 4) {
                                $ins_tp_4bp++;
                        }
			print INSTP $chr, "\t", $begin, "\t", $end, "\tTP\t", $len, "\n";
                }
                else { 
			$ins_fn++; 
			if ($len >= 4) {
				$ins_fn_4bp++
			}
			print INSFN $chr, "\t", $begin, "\t", $end, "\tFN\t", $len, "\n";
		}
        }
}


$snps_tp += $mnvariation->{snp_num} - $mnvariation->{snp_num_priv};
$del_tp += $mnvariation->{del_num} - $mnvariation->{del_num_priv};
$ins_tp += $mnvariation->{ins_num} - $mnvariation->{ins_num_priv};

# print out this stuff
open OUTPUT, ">mn_cmp_results.txt";
print OUTPUT "SNPs:";
print OUTPUT " total:\t", $snps_tp+$snps_fn, "\n";
print OUTPUT " TP:\t", $snps_tp, "\n";
print OUTPUT " FP:\t", $snps_fp, "\n";
print OUTPUT " FN:\t", $snps_fn, "\n";

print OUTPUT "Insertions:";
print OUTPUT " total:\t", $ins_tp+$ins_fn, "\n";
print OUTPUT " TP:\t", $ins_tp, "\n";
print OUTPUT " FP:\t", $ins_fp, "\n";
print OUTPUT " FN:\t", $ins_fn, "\n";

print OUTPUT "Deletions:";
print OUTPUT " total:\t", $del_tp+$del_fn, "\n";
print OUTPUT " TP:\t", $del_tp, "\n";
print OUTPUT " FP:\t", $del_fp, "\n";
print OUTPUT " FN:\t", $del_fn, "\n";

print OUTPUT "Insertions: (>= 4bp)";
print OUTPUT " total:\t", $ins_tp_4bp+$ins_fn_4bp, "\n";
print OUTPUT " TP:\t", $ins_tp_4bp, "\n";
print OUTPUT " FP:\t", $ins_fp_4bp, "\n";
print OUTPUT " FN:\t", $ins_fn_4bp, "\n";

print OUTPUT "Deletions: (>= 4bp)";
print OUTPUT " total:\t", $del_tp_4bp+$del_fn_4bp, "\n";
print OUTPUT " TP:\t", $del_tp_4bp, "\n";
print OUTPUT " FP:\t", $del_fp_4bp, "\n";
print OUTPUT " FN:\t", $del_fn_4bp, "\n";











