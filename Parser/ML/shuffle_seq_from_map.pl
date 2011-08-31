#!/usr/bin/perl
#written by korbinian schneeberger

use strict;
use warnings;
use Getopt::Long;

my %CMD;
my $num;
my $map;
my %ID = ();


GetCom();

open FILE, $map;
my $count_id = 0;

# Read in Read ids
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (!defined($ID{$a[3]})) {
		$count_id++;
	}
	$ID{$a[3]} = 0;
}
close FILE;

print STDERR "got read ids\n";

#  Print new map.list file
open FILE, $map;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (defined($ID{$a[3]}) {
		if ($ID{$a[3]} == 1) {
			print $line;
		}
	}
	else {
		my $rand = rand();
		if ($rand <= ($num / $count_id)) {
                	$num--;
        	        $ID{$a[3]} = 1;
			print $line;
	        }
		else {
			$ID{$a[3]} = 0;
		}
	        $count_id--;
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
	die("Please specify number of entries\n") unless $CMD{num};
  
	$map = $CMD{file};
	$num = $CMD{num};

}
