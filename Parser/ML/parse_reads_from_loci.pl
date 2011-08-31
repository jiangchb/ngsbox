#!/usr/bin/perl
####################################################################################
#Author 	Korbinian Schneeberger 
#Date 		05/15/07
#Version	0.1
#Input		map.list
#Function	Parse out reads from a specified loci
####################################################################################

use strict;
use warnings;
use Getopt::Long;

my %CMD = ();
my $file;
my $begin;
my $end;
my $chr = 1;
my $idx;

GetCom();

open FILE, $file or die "Cannot open $file\n";
open IDX, $idx or die "Cannot open $idx\n";

my $start_pos = $begin;
$start_pos = 1 if $start_pos < 1;


############# Read in index file an jump FILE to right position ############
my $jump = 0;
SETTING: while (my $l = <IDX> ) {
	my @a = split " ", $l;
	if ($a[0] > $chr or ($a[0] == $chr and $a[1] >= $start_pos)) {
		last SETTING;
	} else {
		$jump = $a[2];
	}
}
seek(FILE, $jump, 0);

print STDERR "GOT IT $jump\n";

# set file to region of interest and read in first fragment
INIT: while (<FILE>) {
	my @a = split " ", $_;
	if ($a[0] == $chr and $a[1] >= $start_pos) {
		my $id = $a[3];
		my $seq = parse_seq($a[2]);
		print ">$id\n$seq\n";
		last INIT;
	}
}

print STDERR "FILE SET\n";

# Parse until the end of region of interest

REG: while (<FILE>) {

	# Test for end of region
	my @a = split " ";
	if ($chr != $a[0] or $a[1] > $end) {
		last REG;
	}

	my $id = $a[3];
        my $seq = parse_seq($a[2]);
        print ">$id\n$seq\n";

}


sub parse_seq {
	my ($seq) = @_;

	my $read = "";

	for (my $i=0; $i<length($seq); $i++) {
		if (substr($seq, $i, 1) ne "[") {
			$read .= substr($seq, $i, 1);
		} else {
			$i+=2;
			$read .= substr($seq, $i, 1);
			$i+=1;
		}
	}

	return($read);
}


sub GetCom {

  my @usage = ("\nUsage: $0 --file=<file> --idx=<file> --chr=<int> --begin=<int> --end=<int>

required:
--file\t\tmap.list
--begin\t\tbegin position
--end\t\tend position
--chr\t\tMouse on Mars
--idx\t\tidx file

\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "file=s", "chr=s", "begin=s", "end=s");

	die("Please specify an input file\n") unless $CMD{file};
	die("Please specify begin\n") unless $CMD{begin};
	die("Please specify end\n") unless $CMD{end};
	die("Please specify chr\n") unless $CMD{chr};

	$file = $CMD{file};
	$begin = $CMD{begin};
        $end = $CMD{end};
	$chr = $CMD{chr};

	if (defined($CMD{idx})) {
                $idx = $CMD{idx};
        } else {
                $idx = $file.".idx";
        }	

	return(0);
}

exit(0);
