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
#  Module: Analysis::Assembly::QualityCaller::call_ctg_quality_new.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "\n$0 scaffold_file quality_reference_file quality_variant_file qual_filename masked_filename\n\n";

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

# what variation are considered:
my $min_var_freq = 0.2;
my $min_var_support = 2;

# how to treat ambiguous positions
my $N_extension = 10;
my $N_gap_closing_size = 100;

# quality parameters
my $masking_threshold = 15;
my $min_high_qual_calls = 300;



#####################################
# Statistics

my $passed_scaffolds = 0;
my $dismissed_scaffolds = 0;
my $dismissed500_scaffolds = 0;

my $N_passed_scaffold = 0;
my $N_passed_masked_qual = 0;
my $N_passed_masked_1 = 0;
my $N_passed_masked_2 = 0;
my $N_passed_masked_3 = 0;
my $N_passed_masked_4 = 0;

my $N_dismissed_scaffold = 0;
my $N_dismissed_masked_qual = 0;
my $N_dismissed_masked_1 = 0;
my $N_dismissed_masked_2 = 0;
my $N_dismissed_masked_3 = 0;
my $N_dismissed_masked_4 = 0;

my $N_dismissed500_scaffold = 0;
my $N_dismissed500_masked_qual = 0;
my $N_dismissed500_masked_1 = 0;
my $N_dismissed500_masked_2 = 0;
my $N_dismissed500_masked_3 = 0;
my $N_dismissed500_masked_4 = 0;


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
print STDERR "#:\t", $passed_scaffolds, "\n";
print STDERR "N scaffold:\t\t$N_passed_scaffold\n";
print STDERR "N masked qual:\t\t$N_passed_masked_qual\n";
print STDERR "N masked 1:\t\t$N_passed_masked_1\n";
print STDERR "N masked 2:\t\t$N_passed_masked_2\n";
print STDERR "N masked 3:\t\t$N_passed_masked_3\n";
print STDERR "N masked 4:\t\t$N_passed_masked_4\n";

print STDERR "Scaffolds dismissed:\n";
print STDERR "#:\t", $dismissed_scaffolds, "\n";
print STDERR "N scaffold:\t\t$N_dismissed_scaffold\n";
print STDERR "N masked qual:\t\t$N_dismissed_masked_qual\n";
print STDERR "N masked 1:\t\t$N_dismissed_masked_1\n";
print STDERR "N masked 2:\t\t$N_dismissed_masked_2\n";
print STDERR "N masked 3:\t\t$N_dismissed_masked_3\n";
print STDERR "N masked 4:\t\t$N_dismissed_masked_4\n";

print STDERR "Scaffolds dismissed500:\n";
print STDERR "#:\t", $dismissed500_scaffolds, "\n";
print STDERR "N scaffold:\t\t$N_dismissed500_scaffold\n";
print STDERR "N masked qual:\t\t$N_dismissed500_masked_qual\n";
print STDERR "N masked 1:\t\t$N_dismissed500_masked_1\n";
print STDERR "N masked 2:\t\t$N_dismissed500_masked_2\n";
print STDERR "N masked 3:\t\t$N_dismissed500_masked_3\n";
print STDERR "N masked 4:\t\t$N_dismissed500_masked_4\n";


exit(0);

#####################################
# parse 

sub parse_scaff {
	my ($id, $seq) = @_;	

	# init seq and qual
	my @qual = ();
	my @seq = split "", $seq;

	# set short id
        my @a = split " ", $id;
        my $srt_id = substr($a[0], 1, length($a[0])-1);
        $srt_id =~ s/Scaffold_//;

	# calculate per base quality
	# N = 0
	set_quality($srt_id, \@seq, \@qual);

	# mask scaffolds at contigs of low quality
	# Split = -1
	mask_by_quality_contig(\@seq, \@qual);

	# mask scaffolds at regions of low quality
	# Split = -2
	mask_by_quality_region(\@seq, \@qual);

	# mask short tracks between masked regions
	# Masking = -3
	extend_combine_low_quality_regions(\@seq, \@qual);

	# Print
	# mask by quality missing... ($masking_threshold) (Im print erst?)
	print_seq(\@seq, \@qual);


}

sub print_seq_debug {
	my ($seq_ref, $qual_ref) = @_;

	# Print 
	print QFILE "$id\n";
	print MFILE "$id\n";
	my $pos_c = 0;

	my $print_line_number = 1;
	my $ln = 0;

	for (my $i = 0; $i < @{$seq_ref}; $i++) {

		# line breaks
		if (($pos_c % 50) == 0) {
			if ($pos_c!=0) {
				print QFILE "\n";
				print MFILE "\n";
			}
			if ($print_line_number == 1) {
                                print QFILE $ln*50+1, ":\t";
                                print MFILE $ln*50+1, ":\t";
                                $ln++;
                        }
		}

		# values
		if (${$qual_ref}[$i] < 10 and ${$qual_ref}[$i] >= 0) {
			print QFILE " ";
		}
		print QFILE ${$qual_ref}[$i]." ";
		print MFILE ${$seq_ref}[$i];		

		$pos_c++;
	}
	print QFILE "\n";
	print MFILE "\n";

}

sub print_seq {
        my ($seq_ref, $qual_ref) = @_;

	my $begin = -1;
	my $end = -1;

	# set qual counter and
	# begin
	my $high_qual_counter = 0;
	for (my $i = 0; $i < @{$seq_ref}; $i++) {
		if (${$qual_ref}[$i] >= $masking_threshold) {
                        $high_qual_counter++;
			if ($begin == -1) {
				$begin = $i;
			}
                }
	}

	# set end
	for (my $i = @{$seq_ref}-1; $i>=0; $i--) {
		if (${$qual_ref}[$i] >= $masking_threshold) {
			$end = $i;
			$i = 0; # break loop
		}
	}
		
	# Print 
	if ($high_qual_counter >= $min_high_qual_calls) {
		$passed_scaffolds++;
		print QFILE "$id\n";
		print MFILE "$id\n";
		my $pos_c = 0;
		if ($end > $begin) {
			$N_passed_masked_4 += (@{$seq_ref} - $end) + ($begin - 1);
		}
		else {
			$N_passed_masked_4 += @{$seq_ref}+0;
		}
		for (my $i = $begin; $i <= $end; $i++) {
			if (${$seq_ref}[$i] eq "N") {
				$N_passed_scaffold++;
			}
			if (${$qual_ref}[$i] != -4) {
				# line breaks
				if ($pos_c!=0 and ($pos_c % 40) == 0) {
					print QFILE "\n";
				}
				if ($pos_c!=0 and ($pos_c % 80) == 0) {
					print MFILE "\n";
				}
		
				# print
				# qval
				if (${$qual_ref}[$i]< 10 ) {
					print QFILE " ";
				}
				if (${$qual_ref}[$i]<0) {
					print QFILE "0 ";
				}
				else {
					print QFILE ${$qual_ref}[$i]." ";
				}
				# seq
				if (${$qual_ref}[$i]<$masking_threshold) {
					if (${$seq_ref}[$i] ne "N") {
						if (${$qual_ref}[$i] >= 0) {
							$N_passed_masked_qual++;
						}
						elsif (${$qual_ref}[$i] == -1) {
							$N_passed_masked_1++;
						}
						elsif (${$qual_ref}[$i] == -2) {
                                                        $N_passed_masked_2++;
                                                }
						elsif (${$qual_ref}[$i] == -3) {
                                                        $N_passed_masked_3++;
                                                }
					}
					print MFILE "N";
				}
				else {	
					print MFILE ${$seq_ref}[$i];
				}
		
				$pos_c++;
			}
		}
		print QFILE "\n";
		print MFILE "\n";
	}
	else {
		$dismissed_scaffolds++;
		$dismissed500_scaffolds++ if @{$seq_ref} >= 500;
		print QFILE_DISM "$id\n";
		print MFILE_DISM "$id\n";
		my $pos_c = 0;
                if ($end > $begin) {
                        $N_dismissed_masked_4 += (@{$seq_ref} - $end) + ($begin - 1);
			$N_dismissed500_masked_4 += (@{$seq_ref} - $end) + ($begin - 1) if @{$seq_ref} >= 500;
                }
                else {
                        $N_dismissed_masked_4 += @{$seq_ref}+0;
                        $N_dismissed_masked_4 += @{$seq_ref}+0 if @{$seq_ref} >= 500;
                }
		for (my $i = $begin; $i < $end; $i++) {
			if (${$seq_ref}[$i] eq "N") {
                                $N_dismissed_scaffold++;
                                $N_dismissed500_scaffold++ if @{$seq_ref} >= 500;
                        }
			if (${$qual_ref}[$i] != -4) {
				# lines breaks
				if ($pos_c!=0 and ($pos_c % 40) == 0) {
					print QFILE_DISM "\n";			
				}
		
				if ($pos_c!=0 and ($pos_c % 80) == 0) {
					print MFILE_DISM "\n";
				}
				
				# values
				if (${$qual_ref}[$i] < 10 ) {
					print QFILE_DISM " ";
				}
				if (${$qual_ref}[$i]<0) {
                                        print QFILE_DISM "0 ";
                                }
                                else {
                                        print QFILE_DISM ${$qual_ref}[$i]." ";
                                }
				if (${$qual_ref}[$i]<$masking_threshold) {
					if (${$seq_ref}[$i] ne "N") {
                                                if (${$qual_ref}[$i] >= 0) {
                                                        $N_dismissed_masked_qual++;
							$N_dismissed500_masked_qual++ if @{$seq_ref} >= 500;
                                                }
                                                elsif (${$qual_ref}[$i] == -1) {
                                                        $N_dismissed_masked_1++;
							$N_dismissed500_masked_1++ if @{$seq_ref} >= 500;
                                                }
                                                elsif (${$qual_ref}[$i] == -2) {
                                                        $N_dismissed_masked_2++;
                                                        $N_dismissed500_masked_2++ if @{$seq_ref} >= 500;
                                                }
                                                elsif (${$qual_ref}[$i] == -3) {
                                                        $N_dismissed_masked_3++;
                                                        $N_dismissed500_masked_3++ if @{$seq_ref} >= 500;
                                                }
                                        }
                                        print MFILE_DISM "N";
                                }
                                else {
                                        print MFILE_DISM ${$seq_ref}[$i];
                                }
		
				$pos_c++;
			}
		}
		print QFILE_DISM "\n";
		print MFILE_DISM "\n";
	}


}


sub extend_combine_low_quality_regions {
	my ($seq_ref, $qual_ref) = @_;

	# Extend low quality regions 
	my $n_start = -2;
	my $last_pos = -2;
	for (my $i = 0; $i <= @{$seq_ref}; $i++) {
		# extend N stretches
		if ($i != @{$seq_ref}+0 and (${$qual_ref}[$i] < 15)) {
			$last_pos = $i;
			if ($n_start == -2) {
				$n_start = $i;
			}
		}
		else {
			if ($last_pos != -2 and $last_pos - $n_start + 1>=$N_extension) {
				for (my $j = $n_start-$N_extension; $j < $n_start; $j++) {
					if ($j > 0) {
						if (${$qual_ref}[$j] >= 0) {
							${$qual_ref}[$j] = -3;
							#$n_masked_ext++; # stat
						}
					}
				}
				for (my $j = $last_pos+1; $j < $last_pos+1+$N_extension; $j++) {
					if ($j < @{$seq_ref}) {
						if (${$qual_ref}[$j] >= 0) {
							${$qual_ref}[$j] = -3;
							#$n_masked_ext++; # stat
						}
					}
				}
			}
			$n_start = -2;
			$last_pos = -2;
		}
	}


	# Close gaps between Ns
	my $last_n = -1;
	for (my $i = 0; $i <= @{$seq_ref}; $i++) {
		if ($i==@{$seq_ref} or ${$qual_ref}[$i] < 15) {
			if ($last_n >= $i - $N_gap_closing_size) { 
				for (my $j = $last_n+1; $j < $i; $j++) {
					if (${$qual_ref}[$j] >= 0) {
						${$qual_ref}[$j] = -3;
						#$n_masked_closure++; # stat
					}
				}
			}
			$last_n = $i;
		} 
	}

}

# Search for region in scaffolds that indicate wrong connections.
# This sets the regions that shall be broken to -2 in the quality string
sub mask_by_quality_region {
	my ($seq_ref, $qual_ref) = @_;

	my @pos = ();

	for (my $i = 0; $i < @$seq_ref; $i++) {
		if (${$qual_ref}[$i] < 15 and ${$qual_ref}[$i] >= 0) {
			if (@pos > 0 and $pos[$#pos] < $i - 5) {
				for (my $j = $pos[0]; $j <= $pos[$#pos]; $j++) {
					if (${$qual_ref}[$j] > 0) {
						${$qual_ref}[$j] = -2;
					}
				}
				@pos = ();
			}
			push @pos, $i;
		}
		else {
		}
	}

}

# Search for contigs in scaffolds that should be taken out. The scaffolding
# should be broken(?)
# This sets the regions that shall be broken to -1 in the quality string
sub mask_by_quality_contig {
	my ($seq_ref, $qual_ref) = @_;

	my $count_quality_bases = 0;
	my $count_bases = 0;
	my $start = -1;

	for (my $i = 0; $i < @$seq_ref; $i++) {
		if (${$seq_ref}[$i] eq "N") {
			if ($start != -1) {
				if ($count_quality_bases < 20 or $count_quality_bases/$count_bases <= 0.2) {
					# last region had to little high quality values => mask it (= set quality to -1)
					for (my $j = $start; $j <= $i-1; $j++) {
						if (${$qual_ref}[$j] > 0) {
							${$qual_ref}[$j] = -1;
						}
					}
				}
				$start = -1;
			}
			$count_quality_bases = 0;
			$count_bases = 0;
		}
		else {
			# set the begining of a new contig region after an N stretch
			if ($start == -1) {
				$start = $i;
			}
			# Is it a good or a bad base?
			if (${$qual_ref}[$i] >= 25) {
				$count_quality_bases++;
			}
			$count_bases++;
		}
	}
}

sub set_quality {
	my ($srt_id, $seq_ref, $qual_ref) = @_;

        for (my $i = 0; $i < @$seq_ref; $i++) {
                my $qual;
                my $srt_pos = $i+1;
                $srt_pos = $srt_id."#".$srt_pos;
                # set base quality
                if (${$seq_ref}[$i] eq "N") {
                        ${$qual_ref}[$i] = 0;
                }
                else {
                        # calc with quality caller output
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
                        # set quality
                        ${$qual_ref}[$i] = $qual;
                }
        }

}





# 	# stats
# 	my $n = 0;
# 	my $n_masked_qual = 0;
# 	my $n_masked_ext = 0;
# 	my $n_masked_closure = 0;
# 	my $n_chopped = 0;
# 	my $n_chopped_complete = 0;
# 	my $n_chopped_complete_ctg = 0;
# 
# 	# init output seq
# 	my %q_seq = ();
#         my %m_seq = ();
# 	for (my $i = 0; $i < length($seq); $i++) {
#                 $m_seq{$i} = substr($seq, $i, 1);
# 		if (substr($seq, $i, 1) eq "N") {
# 			$n++; # stat
# 		}
#         }
# 
# 
# 
# 
# 
# 	
# 	# add up stats
# 	if ($high_qual_counter >= $min_high_qual_calls) {
# 		$scaffolds++;
# 		$N += $n;
# 		$N_masked += $n_masked_qual+$n_masked_ext+$n_masked_closure;
# 		$N_masked_qual += $n_masked_qual;
# 		$N_masked_ext += $n_masked_ext;
# 		$N_masked_closure += $n_masked_closure;
# 		$N_chopped += $n_chopped;
#                 $N_chopped_complete += $n_chopped_complete;
#                 $N_chopped_complete_ctg += $n_chopped_complete_ctg;
# 	}
# 	else {
# 		$dismissed_scaffolds++;
# 		$N_dismissed += $n;
# 		$N_masked_dismissed += $n_masked_qual+$n_masked_ext+$n_masked_closure;
# 		$N_masked_qual_dismissed += $n_masked_qual;
# 		$N_masked_ext_dismissed += $n_masked_ext;
# 		$N_masked_closure_dismissed += $n_masked_closure;
# 		$N_chopped_dismissed += $n_chopped;
#                 $N_chopped_complete_dismissed += $n_chopped_complete;
#                 $N_chopped_complete_ctg_dismissed += $n_chopped_complete_ctg;
# 	}
# 
# }

sub max {
	my ($a, $b) = @_;
	return $a if $a >= $b;
	return $b;
}



