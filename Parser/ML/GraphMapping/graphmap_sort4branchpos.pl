#! /usr/bin/perl

use strict;

# Sorting the map.list files at the 23.12 positions correctly:
# This:
# 
# 20.1
# 20.10
# 20.2
# 
# will be:
# 
# 20.1
# 20.2
# 20.10
#



my $usage = "$0 map.list\n";
my $file = shift or die $usage;

open FILE, $file or die "Cannot open file\n";

my $flag = 0;
my %LINE = ();

while (my $line = <FILE>) {
	my ($sample, $chr, $pos) = split " ", $line;
	if ($pos =~ /\./) {
		my ($ref, $strain) = split '\.', $pos;
		$LINE{$strain} .= $line;
		$flag = 1;
	}
	else {
		if ($flag == 1) {
			foreach my $strain (sort{$a <=> $b} keys %LINE) {
				print $LINE{$strain};
			}
			%LINE = ();
		}
		print $line;
		$flag = 0;
	}
}
if ($flag == 1) {
	foreach my $strain (sort{$a <=> $b} keys %LINE) {
        	print $LINE{$strain};
	}
        %LINE = ();
}



