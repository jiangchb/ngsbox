#! /usr/bin/perl
use strict;

my $usage = "$0 length file\n";

my $len = shift or die $usage;
my $file = shift or die $usage;

open FILE, $file or die "Cannot open input file";
my $seq = "";
my $id = "";

while ( <FILE> ) {
	chomp;

	if (substr($_, 0, 1) eq ">") {
		if ($seq ne "") {
			cho($seq);
		}
		$seq = "";
	}
	else {
		$seq .= $_;
	}
}

cho($seq);

sub cho {
	my ($seq) = @_;

	for (my $i = 0; $i < length($seq) - $len; $i++) {
		my $shoreid = 1000000000 + $i;
		print "$shoreid\t" , substr($seq, $i, $len), "\t0\t";
		for(my $j = 0; $j < $len; $j++) {
			print "I";
		}
		print "\n";
	}
}


