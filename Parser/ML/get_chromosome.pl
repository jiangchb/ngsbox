#! /usr/bin/perl
use warnings;
use strict;

my $chr = shift;
my $file = shift;

open FILE, $file or die "cannnot open $file\n";

while(<FILE>) {
	my $line = $_;
	my ($chromosome) = split(/\t/, $line);
	if($chromosome == $chr) {
		print $line;
	}
}

exit(0);
