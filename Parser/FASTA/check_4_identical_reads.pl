#! /usr/bin/perl
use strict;

open FILE, shift;

my %SEQ= ();
my $seq = "";
my $id = "";


### Check fasta file for duplicate sequences
while (<FILE>) {

	chomp();

	if (substr($_, 0, 1) eq ">") {
		my @a = split " ";
		if ($seq ne "") {
			if (defined($SEQ{$seq})) {
				print "Duplicates: ", $SEQ{$seq}, " $id\ņ";
			}
			$SEQ{$seq} .= $id." ";
		}
		$id = substr($a[0], 1, length($a[0])-1);
		$seq = "";
	}
	else {
		$seq.=$_;
	}
}

### Last entry
if ($seq ne "") {
	if (defined($SEQ{$seq})) {
		print "Duplicates: ", $SEQ{$seq}, " $id\ņ";
	}
	$SEQ{$seq} .= $id." ";
}
