#! /usr/bin/perl
use warnings;
use strict;

my $chr_list = shift or die "please specify chromosome list\n";
my $stop = shift;
my $file = shift or die "please specify map.list file\n";

if(! defined $stop) {$stop = -1};

my @chr_tmp = split(",", $chr_list);
my %chr = ();
foreach(@chr_tmp) { $chr{$_} = 1; }

open FILE, $file or die "cannnot open $file\n";

while(<FILE>) {
	my @entries = split(/\t/, $_);
	if( exists $chr{$entries[0]} ) {
		print $_;
	}
	if($entries[0] == $stop) { last; }
}

exit(0);
