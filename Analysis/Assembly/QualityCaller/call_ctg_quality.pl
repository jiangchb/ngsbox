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
#  Module: Analysis::Assembly::QualityCaller::call_ctg_quality.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "\n$0 Q_Threshold Min_HighQual_Calls scaffold_file quality_reference_file quality_variant_file qual_filename masked_filename\n\n";
my $masking_threshold = shift;
my $min_high_qual_calls = shift;
my $s_file = shift or die $usage;
my $r_file = shift or die $usage;
my $v_file = shift or die $usage;
my $q_file = shift or die $usage;
my $m_file = shift or die $usage;

# Input
open SFILE, $s_file or die "Cannot open file ".$s_file."\n";
open RFILE, $r_file or die "Cannot open file ".$r_file."\n";
open VFILE, $v_file or die "Cannot open file ".$v_file."\n";
# Output
open QFILE, ">".$q_file or die "Cannot open file ".$q_file."\n";
open MFILE, ">".$m_file or die "Cannot open file ".$m_file."\n";
open QFILE_DISM, ">".$q_file.".dropped" or die "Cannot open file ".$q_file."\n";
open MFILE_DISM, ">".$m_file.".dropped" or die "Cannot open file ".$m_file."\n";

my %QREF = ();
my %QVAR = ();

#####################################
# Adjustable parameter

my $min_var_freq = 0.25;
my $min_var_support = 2;

my $N_extension = 5;
my $N_extension_region_size = 10;

my $N_gap_closing_size = 10;

#####################################
# Statistics

my $scaffolds = 0;
my $dismissed_scaffolds = 0;

my $N = 0;
my $N_masked = 0;
my $N_masked_qual = 0;
my $N_masked_ext = 0;
my $N_masked_closure = 0;
my $N_chopped = 0;
my $N_chopped_complete = 0;
my $N_chopped_complete_ctg = 0;

my $N_dismissed = 0;
my $N_masked_dismissed = 0;
my $N_masked_qual_dismissed = 0;
my $N_masked_ext_dismissed = 0;
my $N_masked_closure_dismissed = 0;
my $N_chopped_dismissed = 0;
my $N_chopped_complete_dismissed = 0;
my $N_chopped_complete_ctg_dismissed = 0;


#####################################
# Read in quality ref
print STDERR "Reading ref qual...\n";
while (my $l = <RFILE>) {
	my @a = split " ", $l;
	$a[1] =~ s/Scaffold_//;
	$QREF{$a[1]."#".$a[2]} = $a[5]; 
}
close RFILE;
print STDERR "done.\n";

#####################################
# Read in quality var
print STDERR "Reading var qual...\n";
while (my $l = <VFILE>) {
        my @a = split " ", $l;
	if ($a[7] >= $min_var_freq and $a[6] >= $min_var_support) {
	        $a[1] =~ s/Scaffold_//;
        	$QVAR{$a[1]."#".$a[2]} = $a[5];
	}
}
close VFILE;
print STDERR "done.\n";

#####################################
# Parse scaffolds
my $seq = "";
my $id = "";
while (my $l = <SFILE>) {
	chomp($l);
	if (substr($l, 0, 1) eq ">") {
		if ($seq ne "") {
			#print STDERR "Print $id\n";
			parse_scaff($id, $seq);
		}
		$seq = "";
		$id = $l;
	}
	else {
		$seq .= $l;
	}
}
if ($seq ne "") {
	#print STDERR "Print $id\n";
	parse_scaff($id, $seq);
}

#####################################
# Print stats

print STDERR "Scaffolds passed:\n";
print STDERR "#:", $scaffolds, "\n";
print STDERR "N:\t\t\t\t$N\n"; 
print STDERR "N masked:\t\t\t$N_masked\n"; 
print STDERR "N masked qual:\t\t\t$N_masked_qual\n"; 
print STDERR "N masked ext:\t\t\t$N_masked_ext\n"; 
print STDERR "N masked closure:\t\t$N_masked_closure\n"; 
print STDERR "N chopped:\t\t\t$N_chopped\n";
print STDERR "  N chopped complete Ns:\t$N_chopped_complete\n"; 
print STDERR "  N chopped complete num:\t$N_chopped_complete_ctg\n";

print STDERR "Scaffolds dismissed:\n";
print STDERR "#:", $dismissed_scaffolds, "\n";
print STDERR "N:\t\t\t\t$N_dismissed\n";
print STDERR "N masked:\t\t\t$N_masked_dismissed\n";
print STDERR "N masked qual:\t\t\t$N_masked_qual_dismissed\n";
print STDERR "N masked ext:\t\t\t$N_masked_ext_dismissed\n";
print STDERR "N masked closure:\t\t$N_masked_closure_dismissed\n";
print STDERR "N chopped:\t\t\t$N_chopped_dismissed\n";
print STDERR "  N chopped complete Ns:\t$N_chopped_complete_dismissed\n";
print STDERR "  N chopped complete ctg:\t$N_chopped_complete_ctg_dismissed\n";

exit(0);

#####################################
# parse 

sub parse_scaff {
	my ($id, $seq) = @_;	

	# stats
	my $n = 0;
	my $n_masked_qual = 0;
	my $n_masked_ext = 0;
	my $n_masked_closure = 0;
	my $n_chopped = 0;
	my $n_chopped_complete = 0;
	my $n_chopped_complete_ctg = 0;

	# set short id
	my @a = split " ", $id;
	my $srt_id = substr($a[0], 1, length($a[0])-1);
	$srt_id =~ s/Scaffold_//;
	
	# init output seq
	my %q_seq = ();
        my %m_seq = ();
	for (my $i = 0; $i < length($seq); $i++) {
                $m_seq{$i} = substr($seq, $i, 1);
		if (substr($seq, $i, 1) eq "N") {
			$n++; # stat
		}
        }

	# Extend N regions 
	my $n_start = -2;
	my $last_pos = -2;
	for (my $i = 0; $i <= length($seq); $i++) {
		# extend N stretches
		if ($i != length($seq) and substr($seq, $i, 1) eq "N") {
			$last_pos = $i;
			if ($n_start == -2) {
				$n_start = $i;
			}
		}
		else {
			if ($last_pos != -2 and $last_pos - $n_start + 1 >= $N_extension_region_size) {
				for (my $j = $n_start-$N_extension; $j < $n_start; $j++) {
					if ($j > 0) {
						if ($m_seq{$j} ne "N") {
							$m_seq{$j} = "N";
							$n_masked_ext++; # stat
						}
					}
				}
				for (my $j = $last_pos+1; $j < $last_pos+1+$N_extension; $j++) {
					if ($j < length($seq)) {
						if ($m_seq{$j} ne "N") {
							$m_seq{$j} = "N";
							$n_masked_ext++; # stat
						}
					}
				}
			}
			$n_start = -2;
			$last_pos = -2;
		}
	}

	# Mask by quality
	for (my $i = 0; $i < length($seq); $i++) {
		my $qual;
		my $srt_pos = $i+1;
		$srt_pos = $srt_id."#".$srt_pos;
		# set base quality
		if ($m_seq{$i} eq "N") {
			$q_seq{$i} = 0;
		}
		else {
			# set by quality caller output
			if (defined($QREF{$srt_pos}) and not defined($QVAR{$srt_pos})) {
				$qual = $QREF{$srt_pos}
			}
			elsif (defined($QVAR{$srt_pos}) and not defined ($QREF{$srt_pos})) {
				$qual = 0;
			}
			elsif (defined($QREF{$srt_pos}) and defined($QVAR{$srt_pos})) {
				$qual = max(0, $QREF{$srt_pos} - $QVAR{$srt_pos});
			}
			else {
				$qual = 0;
			}
			# print
			if ($qual >= $masking_threshold) {
				$m_seq{$i} = substr($seq, $i, 1);
				$q_seq{$i} = $qual;
			}
			else {
				$m_seq{$i} = "N";
				$q_seq{$i} = 0;
				$n_masked_qual++;
			}
		}
	}


	# Close gaps between Ns
	my $last_n = -1;
	for (my $i = 0; $i <= length($seq); $i++) {
		if ($i==length($seq) or $m_seq{$i} eq "N") {
			if ($last_n >= $i - $N_gap_closing_size) { 
				for (my $j = $last_n+1; $j < $i; $j++) {
					if ($m_seq{$j} ne "N") {
						$m_seq{$j} = "N";
						$q_seq{$j} = 0;
						$n_masked_closure++; # stat
					}
				}
			}
			$last_n = $i;
		} 
	}

	# set beginning and end
	my $begin = 0;
	my $end = length($seq)-1;

	for (my $i = 0; $i < length($seq); $i++) {
		if ($m_seq{$i} eq "N") {
			$begin++;
		}
		else {
			$i = length($seq); # break loop
		}
	}

	for (my $i = length($seq)-1; $i >= 0; $i--) {
                if ($m_seq{$i} eq "N") {
                        $end--;
                }
                else {
                        $i = -1; # break loop
                }
        }
	
	if ($begin < $end) {
		$n_chopped += $begin + (length($seq)-$end);
	}
	else {
		$n_chopped += length($seq);
		$n_chopped_complete += length($seq);
		$n_chopped_complete_ctg++;
	}

	# set qual counter
	my $high_qual_counter = 0;
	for (my $i = $begin; $i < $end; $i++) {
		if ($q_seq{$i} >= $masking_threshold) {
                        $high_qual_counter++;
                }
	}

	# Print 
	if ($high_qual_counter >= $min_high_qual_calls) {
		print QFILE "$id\n";
		print MFILE "$id\n";
		my $pos_c = 0;
		for (my $i = $begin; $i <= $end; $i++) {
			# line breaks
			if ($pos_c!=0 and ($pos_c % 40) == 0) {
				print QFILE "\n";
			}
			if ($pos_c!=0 and ($pos_c % 80) == 0) {
				print MFILE "\n";
			}
	
			# values
			if ($q_seq{$i} < 10 ) {
				print QFILE " ";
			}
			print QFILE $q_seq{$i}." ";
			print MFILE $m_seq{$i};		
	
			$pos_c++;
		}
		print QFILE "\n";
		print MFILE "\n";
	}
	else {
		print QFILE_DISM "$id\n";
		print MFILE_DISM "$id\n";
		my $pos_c = 0;
		for (my $i = $begin; $i <= $end; $i++) {
			# lines breaks
			if ($pos_c!=0 and ($pos_c % 40) == 0) {
				print QFILE_DISM "\n";			
			}
	
			if ($pos_c!=0 and ($pos_c % 80) == 0) {
				print MFILE_DISM "\n";
			}
			
			# values
			if ($q_seq{$i} < 10 ) {
				print QFILE_DISM " ";
			}
			print QFILE_DISM $q_seq{$i}." ";
			print MFILE_DISM $m_seq{$i};
	
			$pos_c++;
		}
		print QFILE_DISM "\n";
		print MFILE_DISM "\n";
	}
	
	# add up stats
	if ($high_qual_counter >= $min_high_qual_calls) {
		$scaffolds++;
		$N += $n;
		$N_masked += $n_masked_qual+$n_masked_ext+$n_masked_closure;
		$N_masked_qual += $n_masked_qual;
		$N_masked_ext += $n_masked_ext;
		$N_masked_closure += $n_masked_closure;
		$N_chopped += $n_chopped;
                $N_chopped_complete += $n_chopped_complete;
                $N_chopped_complete_ctg += $n_chopped_complete_ctg;
	}
	else {
		$dismissed_scaffolds++;
		$N_dismissed += $n;
		$N_masked_dismissed += $n_masked_qual+$n_masked_ext+$n_masked_closure;
		$N_masked_qual_dismissed += $n_masked_qual;
		$N_masked_ext_dismissed += $n_masked_ext;
		$N_masked_closure_dismissed += $n_masked_closure;
		$N_chopped_dismissed += $n_chopped;
                $N_chopped_complete_dismissed += $n_chopped_complete;
                $N_chopped_complete_ctg_dismissed += $n_chopped_complete_ctg;
	}

}

sub max {
	my ($a, $b) = @_;
	return $a if $a >= $b;
	return $b;
}



