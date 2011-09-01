#! /usr/bin/perl
use strict;
use warnings;


my $usage = "\n$0 background,background,... variants\n\n";

my $bg = shift or die $usage;
my $var = shift or die $usage;

my %BG = ();
my %VAR = ();

my @files = split ",", $bg;

foreach my $file (@files) {

	open FILE, $file or die $usage;

	while (my $line = <FILE>) {
		my @a = split " ", $line;
		$BG{$a[1]."#".$a[2]} = 1;
	}

	close FILE;
}


open FILE, $var or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (not defined($BG{$a[1]."#".$a[2]})) {
		my $id = ($a[1] * 100000000) + $a[2];
		$VAR{$id} = $line;
	}

}

foreach my $key (sort {$a <=> $b} keys %VAR) {
	print $VAR{$key};
}



