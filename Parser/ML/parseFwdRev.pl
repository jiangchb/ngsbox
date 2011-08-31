#! /usr/bin/perl
use strict;
use warnings;

my $file = shift;	# map.list file from SHORE or PGSP

open FILE, $file or die "Cannot open infile\n";
open FWD, ">$file.fwd";
open REV, ">$file.rev";

while(<FILE>) {
	my @a = split " ", $_;
	if ($a[4] eq "D") {
		print FWD $_;
	} else {
		print REV $_;
	}
}


