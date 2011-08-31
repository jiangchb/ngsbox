#!/usr/bin/perl

# --------------------------------------------------------------------
#  Written by Stephan Ossowski and Korbinian Schneeberger
#  --------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;
use FindBin;


### Command line parameters
my $consensus;
my $window_size;
my $depth;
my $chr;
my $pool;

my %CMD;
GetCom();


### Parse consensus file
my $chr_size = 0;
my %COV = ();
my %EXP = ();
#my %REP = ();
#my %GC  = ();

open CONS, $consensus or die "Cannot open consensus file\n";

while( <CONS> ) {

	my @a = split("\t", $_);

	if($a[0] eq $chr) {
		$COV{$a[1]} = $a[3];
		$EXP{$a[1]} = $a[47];
		#$GC{$a[1]}  = $a[45];
		#$REP{$a[1]} = $a[10];
		$chr_size = $a[1];
	}
	elsif($a[0] > $chr) {
		last;
	}
}

close CONS or die "Consensus file won't close\n";
print STDERR "Finished reading consensus file\n";

### Create sliding windows --------------------------------------------------------------------
my %DEL = ();
my %CNV = ();

my $sum_cov;
my $sum_exp;
for(my $x = 1; $x <= $window_size; $x++) {
	if(exists $COV{$x}) {
		$sum_cov += $COV{$x};
		$sum_exp += $EXP{$x};
	}
	else {
		$sum_cov += $depth;
		$sum_exp += $depth;
	}
}

for (my $x = 1; $x < ($chr_size - $window_size - 1); $x++) {

	my $ratio = $sum_cov / $sum_exp;

	if( $ratio >= 1.8 ) {
		$CNV{$x} = $ratio;
	}
	elsif( $ratio <= 0.3 ) {
		$DEL{$x} = $ratio;
	}

	$sum_cov = $sum_cov - $COV{$x} + $COV{$x + $window_size};
	$sum_exp = $sum_exp - $EXP{$x} + $EXP{$x + $window_size};
}
print STDERR "Finished sliding window analysis\n";


### Concatenate deletions
my $last_start = -1;
my $last_pos = -1;
my $win_ratio_sum = 0;
my $win_count = 0;

open DELOUT, ">DEL_chr$chr.txt" or die "Cannot open deletion file\n";

foreach my $pos (sort {$a<=>$b} keys %DEL) {

	if( $pos > ($last_pos + $window_size) ) {
		if( $last_pos != -1) {
			 my $end = $last_pos + $window_size;
			 my $length = $end - $last_start + 1;
			 my $win_ratio = $win_ratio_sum / $win_count;
			 print DELOUT "$pool\t$chr\t$last_start\t$end\t$length\t$win_ratio\n";
		}
		$last_start = $pos;
		$last_pos = $pos;
		$win_ratio_sum = $DEL{$pos};
		$win_count = 1;
	}
	else {
		$last_pos = $pos;
		$win_ratio_sum += $DEL{$pos};
		$win_count++;
	}
}
close DELOUT or die;
print STDERR "Finished printing deletions\n";


### Concatenate CNVs
$last_start = -1;
$last_pos = -1;
$win_ratio_sum = 0;
$win_count = 0;

open CNVOUT, ">CNV_chr$chr.txt" or die "Cannot open CNV file\n";

foreach my $pos (sort {$a<=>$b} keys %CNV) {

	if( $pos > ($last_pos + $window_size) ) {
		if( $last_pos != -1) {
			my $end = $last_pos + $window_size;
			my $length = $end - $last_start + 1;
			my $win_ratio = $win_ratio_sum / $win_count;
			print CNVOUT "$pool\t$chr\t$last_start\t$end\t$length\t$win_ratio\n";
		}
		$last_start = $pos;
		$last_pos = $pos;
		$win_ratio_sum = $CNV{$pos};
		$win_count = 1;
	}

	else {
		$last_pos = $pos;
		$win_ratio_sum += $CNV{$pos};
		$win_count++;
	}
}
close CNVOUT or die;
print STDERR "Finished printing CNVs\n";

exit;	

### Min/Max calculation
sub min {
	my ($a, $b) = @_;
	return $a if $a < $b;
	return $b;
}
sub max {
	my ($a, $b) = @_;
	return $a if $a > $b;
	return $b;
}


### Read command line parameters --------------------------------------------------------------
sub GetCom {

	my @usage = ("$0

--consensus      STRING   SHORE consensus file
--depth          DOUBLE   Estimated sequencing depth
--size           INT      Windowsizes
--chrom          INT      Chromosome
--pool           STRING   Chromosome
\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "consensus=s", "depth=f", "size=s", "chrom=s", "pool=s");
				
        die("Please specify consensus file\n") unless defined($CMD{consensus});
        die("Please specify depth\n") unless defined($CMD{depth});
        die("Please specify size\n") unless defined($CMD{size});
	die("Please specify chr\n") unless defined($CMD{chrom});
	die("Please specify pool\n") unless defined($CMD{pool});

        $consensus   = $CMD{consensus};
	$depth       = $CMD{depth};
	$window_size = $CMD{size};
	$chr         = $CMD{chrom};
	$pool        = $CMD{pool};
}

