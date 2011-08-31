#! /usr/bin/perl
use warnings;
use strict;

my $map_list = shift;
my $pe_list = shift;

### Get PE types
my %pe_types = ();
my @tmp = split(",", $pe_list);
foreach ( @tmp ) { $pe_types{$_} = 1; }

### Get alignments with specified PE types
open FILE, $map_list or die "cannnot open $map_list\n";
while (<FILE>) {
	my @a = split("\t", $_);
	if ( exists $pe_types{$a[9]} ) {
		print $_;
	}
}

exit(0);
