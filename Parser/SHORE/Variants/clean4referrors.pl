#! /usr/bin/perl
use strict;

my $usage = "$0 quality_variants.txt referrors.txt\n";

my $qual_file = shift or die $usage;
my $err_file = shift or die $usage;

my %ERR = ();

open FILE, $err_file or die "Cannot open file\n";

while (my $line = <FILE>) {
	my @a = split " ", $line;
	$ERR{$a[0]."#".$a[1]} = 1;
}

close FILE;

open FILE, $qual_file or die "Cannot open file\n";
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (not defined($ERR{$a[1]."#".$a[2]})) {
		print $line;
	}
}


