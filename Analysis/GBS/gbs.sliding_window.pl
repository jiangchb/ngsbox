#! /usr/bin/perl

## outout file:

#WINDOW_WISE:
#CHR	WINSTART	WINEND	CALL	COUNT_P1	COUNT_P2
#
#BREAKPOINT_WISE:
#CHR	START	END	CALL


#DOWNSTREAM ANALYSIS:
#* OUTPUT RECOMBINATION BREAK NUMS

use strict;

my $usage = "$0 consensus filtered_marker chrsizes\n";
my $f_cons = shift or die $usage;
my $f_marker = shift or die $usage;
my $f_chrsizes = shift or die $usage;
print_log();

my %MARKER_A1 = ();
my %MARKER_A2 = ();
my %MARKER_C1 = ();
my %MARKER_C2 = ();
my %ALLELE = ("A" => 0, "C" => 1, "G" => 2, "T" => 3);
my $REPEAT_FILTERED_BASE_CALLS = 1;

my %SLID_WIN_C1 = ();
my %SLID_WIN_C2 = ();
my %SLID_WIN_M = ();
my %SLID_WIN_CALL = ();

my %CHR_SIZE = ();
my $SLID_WIN_SIZE = 100000;
my $SLID_WIN_MIN_MARKER = 10;
my $SLID_WIN_MIN_READS = 15;

my $SLID_WIN_HOM_CUTOFF = 0.9;
my $SLID_WIN_HET_CUTOFF = 0.10;

### Read in markers and allele counts

read_marker($f_marker);
read_counts($f_cons);
read_chrsizes($f_chrsizes);

### Perform sliding window (very simplistic)

foreach my $chr (sort {$a <=> $b} keys %MARKER_C1) {
	foreach my $pos (sort {$a <=> $b} keys %{$MARKER_C1{$chr}}) {
		my $win = int(($pos-1) / $SLID_WIN_SIZE) * $SLID_WIN_SIZE;
		$SLID_WIN_C1{$chr}{$win} += $MARKER_C1{$chr}{$pos};
		$SLID_WIN_C2{$chr}{$win} += $MARKER_C2{$chr}{$pos};	
		$SLID_WIN_M{$chr}{$win}++; 
	}
}

### Perform genotyping

open OUT, "> GBS.sliding_window_genotyping.txt" or die "cannot open output file\n";

foreach my $chr (sort {$a <=> $b} keys %SLID_WIN_C1) {
	for (my $pos = 0; $pos <= $CHR_SIZE{$chr}; $pos += $SLID_WIN_SIZE) {
		
		my $end = min(($pos+$SLID_WIN_SIZE), $CHR_SIZE{$chr});
		my $c1 = 0;
		my $c2 = 0;
		my $call = "N";
		my $num_marker = (defined($SLID_WIN_M{$chr}{$pos})? $SLID_WIN_M{$chr}{$pos} : 0);
	
		if (defined($SLID_WIN_C1{$chr}{$pos})) {
			$c1 = $SLID_WIN_C1{$chr}{$pos};
			$c2 = $SLID_WIN_C2{$chr}{$pos};
			my $cs = $c1 + $c2;
	
			if ($cs < $SLID_WIN_MIN_READS) {
				$call = "N";
			}
			elsif ($num_marker < $SLID_WIN_MIN_MARKER) {
				$call = "N";
			}
			elsif ($c1 >= $SLID_WIN_HOM_CUTOFF*$cs) {
				$call = "1";
			}
			elsif ($c2 >= $SLID_WIN_HOM_CUTOFF*$cs) {
				$call = "2";
			}
			elsif (($c1/$cs) >= 0.5 - $SLID_WIN_HET_CUTOFF and ($c1/$cs) <= 0.5 + $SLID_WIN_HET_CUTOFF) {
				$call = "3"
			}
		}

		print OUT "$chr\t",($pos+1),"\t",$end,"\t$call\t$c1\t$c2\t",($c1+$c2),"\t$num_marker\n";
		$SLID_WIN_CALL{$chr}{$pos} = $call;
	}
}

close OUT;

### Call breakpoints

open OUT, "> GBS.sliding_window_breakpoints.txt";

my %count_recomb_break = ();

foreach my $chr (sort {$a <=> $b} keys %SLID_WIN_CALL) {
	my $start = 0;
	my $pos = 0;
	my $end = 0;
	my $genotype = 0;

	foreach $pos (sort {$a <=> $b} keys %{$SLID_WIN_CALL{$chr}}) {
		my $call = $SLID_WIN_CALL{$chr}{$pos};
		if ($call ne "N") {
			if ($genotype == 0) {
				$genotype = $call;
				$end = $pos + $SLID_WIN_SIZE;
			}
			elsif ($genotype == $call) {
				$end = $pos + $SLID_WIN_SIZE;
			}
			else {
				## report last genotype
				print OUT $chr,"\t",($start+1),"\t",$end,"\t",$genotype,"\n";
				if ($end != $pos) {
					print OUT $chr,"\t",($end+1),"\t",$pos,"\tN\n";
				}
				## set new background
				$genotype = $call;
				$start = $pos;
				$end = $pos + $SLID_WIN_SIZE;

				$count_recomb_break{$chr}++;
			}
		}
	}

	$end = $CHR_SIZE{$chr};

	if ($genotype == 0) { # no genotype was detected throughout chromosome
		print OUT $chr,"\t",$start,"\t",$end,"\tN\n";
	}
	else {
		print OUT $chr,"\t",$start,"\t",$end,"\t",$genotype,"\n";
	} 
	
}

### Output breakpoint summary
open OUT, "> GBS.sliding_window_recomb_num.txt" or die "cannot open log";
my $f = 0;
#foreach my $chr (sort {$a <=> $b} keys %count_recomb_break) {
#	print OUT "\t" if $f != 0;
#	$f = 1;
#	print OUT $chr;
#}
#print OUT "\tall\n";
print OUT "1\t2\t3\t4\t5\tall\n";
$f = 0;
my $sum=0;
#foreach my $chr (sort {$a <=> $b} keys %count_recomb_break) {
for (my $chr=1; $chr<=5; $chr++) {
        print OUT "\t" if $f != 0;
        $f = 1;
	if (defined($count_recomb_break{$chr})) {
	        print OUT $count_recomb_break{$chr};
		$sum+=$count_recomb_break{$chr};
	}
	else {
	        print OUT "0";
	}
}       
print OUT "\t$sum\n";


##################################################
### Helpers

sub min {
	my ($a, $b) = @_;
	return $a if ($a < $b);
	return $b;
}

sub print_log {
	open LOG, "> GBS.sliding_window.log" or die "cannot open log";;
	print LOG "# $0";
	foreach my $a (@ARGV) {
		print LOG " ", $a;
	}
	print LOG "\n";
	close LOG;
}

sub read_chrsizes {
        my ($file) = @_;
        open FILE, $file or die "cannot open $file\n";
        while (<FILE>) {
                my @a = split " ";
                my $chr = $a[0];
                my $size = $a[1];
		$CHR_SIZE{$chr} = $size;
	}
}

sub read_counts {
	my ($file) = @_;
        open FILE, $file or die "cannot open $file\n";
	my $count = 0;
        while (<FILE>) {
		my @a = split " ";
                my $chr = $a[0];
                my $pos = $a[1];

		if (defined($MARKER_A1{$chr}{$pos})) {
			my $base;
	                if ($REPEAT_FILTERED_BASE_CALLS == 0) {
        	                $base = 5-1;
                	}
	                else {
        	                $base = 27-1;
                	}

                	$MARKER_C1{$chr}{$pos} = $a[$base + $ALLELE{$MARKER_A1{$chr}{$pos}}];
                	$MARKER_C2{$chr}{$pos} = $a[$base + $ALLELE{$MARKER_A2{$chr}{$pos}}];
			$count++;
		}
	}
	print STDERR "Counted alleles at markers: $count\n";
}

sub read_marker {
	my ($file) = @_;
	open FILE, $file or die "cannot open $file\n";
	my $count = 0;
	while (<FILE>) {
		my @a = split " ";
		my $chr = $a[1]; 
		my $pos = $a[2];

		$MARKER_A1{$chr}{$pos} = $a[3];
		$MARKER_A2{$chr}{$pos} = $a[4];
		
		$count++;
	}
	print STDERR "Found markers: $count\n";
}
 
