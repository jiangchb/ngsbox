#! /usr/bin/perl

use strict;

my $usage = "reverse_fl.pl fl-file";
my $file = shift or die $usage;


open FILE, $file or die $usage;

while (<FILE>) {
	my @a = split " ", $_;
	print $a[0], "\t";
	print rev_comp($a[1]), "\t";
	for (my $i = 2; $i < @a; $i++) {
		my $s = reverse $a[$i];
		print $s;
		print "\t" if $i != @a -1;
	}
	print "\n";
}



sub rev_comp {
	my $seq = shift;
	my $rev = "";
	
	for (my $i = 0; $i<length($seq); $i++) {
		$rev .= rev_comp_base(substr($seq, $i, 1));
	}
	
	my $s = reverse $rev;

	return $s;
}

sub rev_comp_base {
        my ($base) = @_;

        $base =~ tr/ACTGactgMRVHmrvhKYBDkybd/TGACtgacKYBDkybdMRVHmrvh/;

        return $base;
}




