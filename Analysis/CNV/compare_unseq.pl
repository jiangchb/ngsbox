#!/usr/bin/perl

# --------------------------------------------------------------------
#  Written by Stephan Ossowski and Korbinian Schneeberger
#  --------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;
use FindBin;


### Command line parameters
my $input;
my $control1;
my $control2;
my $chr;

my %CMD;
GetCom();


### Parse input
my %Input  = ();
my %Input_line = ();
open INPUT, $input or die "Cannot open input file\n";
while( <INPUT> ) {
	my @a = split("\t", $_);
	if( ($a[1] eq $chr) && ($a[4] >= 25) ) {
		$Input{$a[2]} = $a[3];
		$Input_line{$a[2]} = $_;
	}
}
close INPUT or die;


### Parse control 1
my %Control1  = ();
open CONTROL1, $control1 or die "Cannot open control1 file\n";
while( <CONTROL1> ) {
	my @a = split("\t", $_);
	if( ($a[1] eq $chr) && ($a[4] >= 5) ) {
		$Control1{$a[2]} = $a[3];
	}
}
close CONTROL1 or die;


### Parse control 1
my %Control2  = ();
open CONTROL2, $control2 or die "Cannot open input file\n";
while( <CONTROL2> ) {
	my @a = split("\t", $_);
	if( ($a[1] eq $chr) && ($a[4] >= 5) ) {
		$Control2{$a[2]} = $a[3];
	}
}
close CONTROL2 or die;


### Detect overlaps
foreach my $ibeg (sort {$a<=>$b} keys %Input) {
	
	my $iend = $Input{$ibeg};

	# check control 1
	my $total_overlap_1 = 0;
	foreach my $cbeg (sort {$a<=>$b} keys %Control1) {
		my $cend = $Control1{$cbeg};
		
		if( ($cbeg >= $ibeg) && ($cbeg <= $iend) ) {
			my $overlap = min( ($iend - $cbeg + 1), ($cend - $cbeg + 1) );
			$total_overlap_1 += $overlap;
		}
		elsif( ($cend >= $ibeg) && ($cend <= $iend) ) {
			$total_overlap_1 += ($cend - $ibeg + 1);
		}
		elsif( ($ibeg >= $cbeg) && ($iend <= $cend) ) {
			$total_overlap_1 += ($iend - $ibeg);
			last;
		}
	}

	# check control 2
	my $total_overlap_2 = 0;
	foreach my $cbeg (sort {$a<=>$b} keys %Control2) {
		my $cend = $Control2{$cbeg};
		
		if( ($cbeg >= $ibeg) && ($cbeg <= $iend) ) {
			my $overlap = min( ($iend - $cbeg + 1), ($cend - $cbeg + 1) );
			$total_overlap_2 += $overlap;
		}
		elsif( ($cend >= $ibeg) && ($cend <= $iend) ) {
			$total_overlap_2 += ($cend - $ibeg + 1);
		}
		elsif( ($ibeg >= $cbeg) && ($iend <= $cend) ) {
			$total_overlap_2 += ($iend - $ibeg + 1);
			last;
		}
	}
	
	# Print if unseq region is not found in at least one control
	if(
		( $total_overlap_1 < ( ($iend - $ibeg + 1) / 5 ) ) ||
		( $total_overlap_2 < ( ($iend - $ibeg + 1) / 5 ) )
	) {
		print $Input_line{$ibeg};
	}
}

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

--input          STRING   SHORE CNV or unseq file
--control1       STRING   SHORE CNV or unseq file
--control2       STRING   SHORE CNV or unseq file
--chrom          INT      Chromosome
\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "input=s", "control1=s", "control2=s", "chrom=s");
				
        die("Please specify input file\n") unless defined($CMD{input});
        die("Please specify control1 file\n") unless defined($CMD{control1});
        die("Please specify control2 file\n") unless defined($CMD{control2});
	die("Please specify chr\n") unless defined($CMD{chrom});

        $input    = $CMD{input};
	$control1 = $CMD{control1};
	$control2 = $CMD{control2};
	$chr      = $CMD{chrom};
}

