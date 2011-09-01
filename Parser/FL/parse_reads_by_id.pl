#! /usr/bin/perl

use strict;
my $usage = "$0 idfile readsfile\n";
my $idfile = shift or die $usage;
my $readfile = shift or die $usage;

open IDFILE, $idfile or die $usage;
open READFILE, $readfile or die $usage;

my %IDS = ();

while (my $line = <IDFILE>) {
	chomp($line);
	$IDS{$line} = 1;
}

close IDFILE;

while (my $line = <READFILE>) {
	my @a = split " ", $line;
	if (defined($IDS{$a[0]})) {
		print $line;
	}
}

close READFILE;

