#!/usr/bin/perl
#written by korbinian schneeberger

use strict;
use warnings;
use Getopt::Long;

my %CMD;
my $prob;
my $map;
my %ID = ();


GetCom();


#  Print new map.list file
open FILE, $map;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (defined($ID{$a[3]})) {
		if ($ID{$a[3]} == 1) {
			print $line;
		}
	}
	else {
		my $rand = rand();
		if ($rand <= $prob) {
        	        $ID{$a[3]} = 1;
			print $line;
	        }
		else {
			$ID{$a[3]} = 0;
		}
	}
}


sub GetCom{

  my @usage = ("Usage: $0 --file=<map.list> --prob=<double> 

required:
--file\t\tmap.list file to be parsed
--prob\t\tProbability that a read is selected for subsampling
\n"); 
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "file=s","prob=s");

	die("Please specify a fl file\n") unless $CMD{file};
	die("Please specify number of entries\n") unless $CMD{prob};
  
	$map = $CMD{file};
	$prob = $CMD{prob};

}
