#! /usr/bin/perl
use warnings;
use strict;

my $file = shift;
open FILE, $file or die "cannnot open $file\n";

while(<FILE>) {
	my $line = $_;
	my @entries = split(/\t/, $line);
	$entries[1]++;

	if( ($entries[0] =~ /\d+/) && ($entries[1] =~ /\d+/) ) {
		print join("\t", @entries);
	}
}

exit(0);
