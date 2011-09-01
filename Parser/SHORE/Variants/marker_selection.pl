#!/usr/bin/perl

use strict;
use warnings;

## adjustments

my $CVP_EXT = 50;

my $MIN_QUA_VAR = 25;
my $MIN_CON_VAR = 0.9;
my $MAX_REP_VAR = 1.1;

my $MIN_QUA_REF = 25;
my $MIN_CON_REF = 0.9;
my $MAX_REP_REF = 1.1;

my $MIN_DISTANCE = 100;

my $COV_PRC = 0.025;

## Usage

my $usage = "$0 
p1_name p1_variants p1_reference p1_CNV p1_copy_variable_position
p2_name p2_variants p2_reference p2_CNV p2_copy_variable_position
\n";

my $p1_name = shift or die $usage;
my $p1_var_file = shift or die $usage;
my $p1_ref_file = shift or die $usage;
my $p1_cnv_file = shift or die $usage;
my $p1_cvp_file = shift or die $usage;

my $p2_name = shift or die $usage;
my $p2_var_file = shift or die $usage;
my $p2_ref_file = shift or die $usage;
my $p2_cnv_file = shift or die $usage;
my $p2_cvp_file = shift or die $usage;

my $out_file = "HQ_marker.".$p1_name."_x_".$p2_name.".txt";
my $out_id_file = "HQ_marker.".$p1_name."_x_".$p2_name.".shore_id.txt";
my $cov1_file = "HQ_cov1.".$p1_name."_x_".$p2_name.".txt";
my $cov2_file = "HQ_cov2.".$p1_name."_x_".$p2_name.".txt";

## log file
my $log_file = "HQ_marker.".$p1_name."_x_".$p2_name.".log";
if (-e $log_file) {
        print STDERR "$log_file already exists. exit.\n";
        exit(1);
}
open LOG, ">$log_file" or die "cannot open marker_selection.log\n";
open COV1F, ">$cov1_file" or die "cannot open $cov1_file\n";
open COV2F, ">$cov2_file" or die "cannot open $cov2_file\n";
print LOG "time perl $0 $p1_name $p1_var_file $p1_ref_file $p1_cnv_file $p1_cvp_file $p2_name $p2_var_file $p2_ref_file $p2_cnv_file $p2_cvp_file\n";

## read in filter regions

my %MASK = ();
my $TOTAL_BASES_MASKED = 0;

read_cnv($p1_cnv_file);
print "Read cnv1\n";
read_cnv($p2_cnv_file);
print "Read cnv2\n";

read_cvp($p1_cvp_file);
print "Read cvp1\n";
read_cvp($p2_cvp_file);
print "Read cvp2\n";

print LOG "# Total masked positions: $TOTAL_BASES_MASKED\n";

## read in marker

my %MARKER1_RAW = ();
my %MARKER1 = ();
my %MARKER2_RAW = ();
my %MARKER2 = ();

read_var(\%MARKER1, \%MARKER1_RAW, $p1_var_file);
print "Read var1\n";
read_var(\%MARKER2, \%MARKER2_RAW, $p2_var_file);
print "Read var2\n";

filter_overlapping_marker(\%MARKER1, \%MARKER2_RAW);
filter_overlapping_marker(\%MARKER2, \%MARKER1_RAW);

## read in reference calls

my %HQ_MARKER = ();

read_ref($p1_ref_file, \%MARKER2, \%HQ_MARKER, 1);
read_ref($p2_ref_file, \%MARKER1, \%HQ_MARKER, 2);

filter_hq_marker(\%HQ_MARKER);

print_hq_marker(\%HQ_MARKER);


#######################################
## subroutines

sub print_hq_marker {
	my ($markerref) = @_;

	open OUT, "> $out_file" or die "cannot open $out_file\n";
	open OUTID, "> $out_id_file" or die "cannot open $out_file\n";

	my $c_chr = 0;

        foreach my $chr (sort{$a cmp $b}keys %{$markerref}) {
		$c_chr++;
                foreach my $pos (sort{$a<=>$b}keys %{${$markerref}{$chr}}) {
			my $val_str = ${${$markerref}{$chr}}{$pos};
                        my @val = split "#", $val_str;
			print OUT $p1_name."_x_".$p2_name."\t$chr\t$pos\t", $val[1], "\t", $val[3], "\t", $val[4], "\n";
			print OUTID $p1_name."_x_".$p2_name."\t$c_chr\t$pos\t", $val[1], "\t", $val[3], "\t", $val[4], "\n";
		}
	}

	close OUT;
	close OUTID;
	
}

sub filter_hq_marker {
	my ($markerref) = @_;

	my @COV1 = ();
	my @COV2 = ();

	foreach my $chr (sort{$a cmp $b}keys %{$markerref}) {
		foreach my $pos (sort{$a<=>$b}keys %{${$markerref}{$chr}}) {
			my $val_str = ${${$markerref}{$chr}}{$pos};
			my @val = split "#", $val_str;
			push @COV1, $val[0];
			push @COV2, $val[2];
		}
	}

	@COV1 = sort {$a<=>$b} @COV1;
	@COV2 = sort {$a<=>$b} @COV2;

	print COV1F join "\n", @COV1;
	print COV2F join "\n", @COV2;

	my ($lower1, $upper1) = get_cutoffs(\@COV1);
	my ($lower2, $upper2) = get_cutoffs(\@COV2);

	print LOG "# coverage1: $lower1 - $upper1\n";
	print LOG "# coverage2: $lower2 - $upper2\n";

	my $c_del = 0;
	my $c_cov_del1 = 0;
	my $c_cov_del2 = 0;
	my $c_dist_del = 0;

	foreach my $chr (sort{$a cmp $b}keys %{$markerref}) {
		my $last_pos = (-1)*$MIN_DISTANCE;
                foreach my $pos (sort{$a<=>$b}keys %{${$markerref}{$chr}}) {
                        my $val_str = ${${$markerref}{$chr}}{$pos};
       	                my @val = split "#", $val_str;
			if ($val[0] <= $lower1 or $val[0] >= $upper1) {
				$c_cov_del1++;
			}
			if ($val[2] <= $lower2 or $val[2] >= $upper2) {
				$c_cov_del2++;
			}
			if ($last_pos + $MIN_DISTANCE >= $pos) {
				$c_dist_del++;
			}
			if ($last_pos + $MIN_DISTANCE >= $pos or $val[0] <= $lower1 or $val[0] >= $upper1 or $val[2] <= $lower2 or $val[2] >= $upper2) {
				delete ${${$markerref}{$chr}}{$pos};
				$c_del++;
			}
			else {
				$last_pos = $pos;
			}
                }
        }

	
	print LOG "#Markers coverage1 out-of-bounds: $c_cov_del1\n";
	print LOG "#Markers coverage2 out-of-bounds: $c_cov_del2\n";
	print LOG "#Markers distance violation: $c_dist_del\n";
	print LOG "#Markers deleted: $c_del\n";

}

sub get_cutoffs {
	my ($arrref) = @_;
	
	my $mem_num = @{$arrref}+0;
	my $prc_num = $mem_num * $COV_PRC;

	my $lower = ${$arrref}[$prc_num];
	my $upper = ${$arrref}[max(0,$mem_num-$prc_num-1)];

	return ($lower, $upper);
}


sub read_ref {
	my ($file, $markerref, $hq_markerref, $parent) = @_;
	open FILE, $file or die "cannot open $file\n";

	my $C_FOUND = 0;
        my $C_MASK = 0;
        my $C_MIN_QUA_REF = 0;
        my $C_MAX_REP_REF = 0;
        my $C_MIN_CON_REF = 0;
        my $C_SEQ_REF = 0;
        my $C_MAR = 0;

        while (<FILE>) {
                my @a = split " ";
                my $chr = $a[1];
                my $pos = $a[2];
		if (defined(${$markerref}{$chr}{$pos})) {
			$C_FOUND++;
	                if (not defined($MASK{$chr}{$pos})) {		
	                        my $qua = $a[5];
        	                my $ref = $a[3];
                	        my $snp = $a[4];
                        	my $rep = $a[8];
	                        my $con = $a[7];
        	                my $sup = $a[6];
                	        my $cov = int($sup / $con);
	                        if ($qua >= $MIN_QUA_REF) {
        	                        if ($rep <= $MAX_REP_REF) {
                	                        if ($con >= $MIN_CON_REF) {
							my @second_allele = split "#", ${$markerref}{$chr}{$pos};
							my $second_allele_snp = $second_allele[1];
                        	                        if ($ref ne "N" and $snp ne $second_allele_snp and ($snp eq "A" or $snp eq "C" or $snp eq "G" or $snp eq "T")) {
								$C_MAR++;
								if ($parent == 1) {
	                                	                        ${$hq_markerref}{$chr}{$pos} = $cov."#".$snp."#".${$markerref}{$chr}{$pos}."#1";
								}
								else {
	                                	                        ${$hq_markerref}{$chr}{$pos} = ${$markerref}{$chr}{$pos}."#".$cov."#".$snp."#2";
								}
								delete ${$markerref}{$chr}{$pos};
                                        	        }
                                                	else {
                                                        	$C_SEQ_REF++;
	                                                }
        	                                }
                	                        else {
                        	                        $C_MIN_CON_REF++;
                                	        }
	                                }
        	                        else {
                	                        $C_MAX_REP_REF++;
                        	        }
	                        }
        	                else {
                	                $C_MIN_QUA_REF++;
                        	}
	                }
	                else {
        	                $C_MASK++;
                	}
		}
        }
        print LOG "#Reading reference allele call from $file\n";
	print LOG "Found marker:\t$C_FOUND\n";
        print LOG "#Filtered:\n";
        print LOG "Region masking:\t$C_MASK\n";
        print LOG "Quality:\t$C_MIN_QUA_REF\n";
        print LOG "Repetivity:\t$C_MAX_REP_REF\n";
        print LOG "Concordance:\t$C_MIN_CON_REF\n";
        print LOG "Sequence:\t$C_SEQ_REF\n";
        print LOG "#Accepted:\t$C_MAR\n";
	
	open TMP, ">".$parent.".not_used_marker" or die "cannot open file";
	foreach my $chr (sort{$a cmp $b}keys %{$markerref}) {
                foreach my $pos (sort{$a<=>$b}keys %{${$markerref}{$chr}}) {
			print TMP $chr, "\t", $pos, "\n";
		}
	}
}

sub filter_overlapping_marker {
	my ($markerref1, $rawmarkerref2) = @_;

	my $num_d = 0;
	my $num_r = 0;
	
	foreach my $chr (sort{$a cmp $b}keys %{$markerref1}) {
                foreach my $pos (sort{$a<=>$b}keys %{${$markerref1}{$chr}}) {
			if (defined(${$rawmarkerref2}{$chr}{$pos})) {
                        	$num_d++;
				delete ${$markerref1}{$chr}{$pos};
			}
			else {
				$num_r++;
			}
                }
        }

	print LOG "#Overlapping HQ marker: $num_d (deleted)\n";
	print LOG "#Marker remaining: $num_r\n";
}


sub read_var {
	my ($markerref, $rawmarkerref, $file) = @_;
        open FILE, $file or die "cannot open $file\n";

        my $C_MASK = 0;
        my $C_MIN_QUA_VAR = 0;
        my $C_MAX_REP_VAR = 0;
        my $C_MIN_CON_VAR = 0;
        my $C_SEQ_VAR = 0;
        my $C_MAR = 0;
	my $C_ALL = 0;

        while (<FILE>) {
                my @a = split " ";
		my $chr = $a[1];
		my $pos = $a[2];
		my $qua = $a[5];
                my $ref = $a[3];
                my $snp = $a[4];
                my $rep = $a[8];
                my $con = $a[7];
                my $sup = $a[6];
                my $cov = int($sup / $con);
		if (not defined($MASK{$chr}{$pos})) {
			if ($qua >= $MIN_QUA_VAR) {
				if ($rep <= $MAX_REP_VAR) {
					if ($con >= $MIN_CON_VAR) {
						if ($ref ne "N" and $snp ne "-") {
							${$markerref}{$chr}{$pos} = $cov."#".$snp;
							$C_MAR++;
						}
						else {
							$C_SEQ_VAR++;
						}
					}
					else {
						$C_MIN_CON_VAR++;
					}
				}
				else {
					$C_MAX_REP_VAR++;
				}
			}
			else {
				$C_MIN_QUA_VAR++;
			}
		}
		else {
			$C_MASK++;
		}
		${$rawmarkerref}{$chr}{$pos} = $cov."#".$snp;
		$C_ALL++;
	}
	print LOG "#Reading markers from $file\n";
	print LOG "All markers:\t$C_ALL\n";
	print LOG "#Filtered:\n";
	print LOG "Region masking:\t$C_MASK\n";
	print LOG "Quality:\t$C_MIN_QUA_VAR\n";
	print LOG "Repetivity:\t$C_MAX_REP_VAR\n";
	print LOG "Concordance:\t$C_MIN_CON_VAR\n";
	print LOG "Sequence:\t$C_SEQ_VAR\n";
	print LOG "#Accepted:\t$C_MAR\n";
}

sub read_cvp {
        my ($file) = @_;
        open FILE, $file or die "cannot open $file\n";
	my $C_MASK = 0;
	my $last_chr = "";
	my $last_pos = -1;
        while (<FILE>) {
                my @a = split " ";
		if ($a[1] ne $last_chr) {
			$last_pos = -1;
		}
                for (my $i = max(1, max($last_pos+1, $a[2]-$CVP_EXT)); $i <= $a[2]+$CVP_EXT; $i++) {
			if (not defined($MASK{$a[1]}{$i})) {
                        	$MASK{$a[1]}{$i} = 1;
                                $C_MASK++;
			}
                }
		$last_chr = $a[1];
		$last_pos = $a[2]+$CVP_EXT;
        }
        close FILE;
	print LOG "#Reading CVPs from $file\n";
        print LOG "Masked:\t$C_MASK\n";
	$TOTAL_BASES_MASKED += $C_MASK;
}


sub read_cnv {
	my ($file) = @_;
	open FILE, $file or die "cannot open $file\n";
	my $C_MASK = 0;
	while (<FILE>) {
		my @a = split " ";
		for (my $i = $a[2]; $i <= $a[3]; $i++) {
			if (not defined($MASK{$a[1]}{$i})) {
				$MASK{$a[1]}{$i} = 1;
				$C_MASK++;
			}
		}
	}
	close FILE; 
	print LOG "#Reading CNVs from $file\n";
        print LOG "Masked:\t$C_MASK\n";
	$TOTAL_BASES_MASKED += $C_MASK;
}


sub max {
	my ($a, $b) = @_;
	return $a if ($a > $b);
	return $b;
}

