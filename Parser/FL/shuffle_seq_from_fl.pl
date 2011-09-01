#!/usr/bin/perl
#written by korbinian schneeberger

use strict;
use warnings;
use Getopt::Long;

my %CMD;
my $num;
my $seq;

GetCom();

while(<FILE>) {
	my $rand = rand();
	if ($rand <= ($num / $seq)) {
      		$num--;
      		print $_;
	} 
	$seq--;
}


sub GetCom{

  my @usage = ("Usage: shuffle_fl.pl --file=<fl file> --num=<number of entries> 

required:
--file\t\tFl file to be parsed
--num\t\tNumber of sequences randomly taken out of the fl file
--perc\t\tPercent of sequences randonly taken out of the fl file
\n"); 
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "file=s","num=s","perc=s");

	die("Please specify a fl file\n") unless $CMD{file};
	die("Please specify number of entries\n") unless ($CMD{num} || $CMD{perc});

	open FILE, $CMD{file} or die "Cannot open file\n";	
	
	my $out = `wc -l $CMD{file}`;
	my @a = split " ", $out;
	$seq = $a[0];

	if($CMD{num}) {
		$num = $CMD{num};
	}
	elsif($CMD{perc}) {
		$num = int($seq * $CMD{perc} / 100);
	}
}
