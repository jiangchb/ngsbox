#! /usr/bin/perl
use strict;

my $usage = "\n$0 ids_file map_list\n\n";
my $readfile = shift or die $usage;
my $file = shift or die $usage;

my %IDS = ();

open FILE, $readfile or die "Cannot open file\n";
while (my $line = <FILE>) {
	chomp($line);
	$IDS{$line} = 1;
}
close FILE;

open FILE, $file or die "Cannot open file\n";
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (defined($IDS{$a[3]})) {
		print $line;
	}
}
