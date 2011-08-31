#! /usr/bin/perl

use strict;
use warnings;


open SNP, $ARGV[0];
open QUAL, $ARGV[1];

my $chr = 0;
my $pos = 1;
my %POS2QUAL = ();

QUALITY: while (my $line = <QUAL>) {
	chomp($line);
	if (substr($line, 0, 1) eq ">") {
		# debug
		if ($line ne ">scaffold_1") {
			last QUALITY;
		}
		print STDERR $line, "\n";
		$pos = 1;	
		$chr++;
	}
	else {
		my @a = split " ", $line;
		for (my $i = 0; $i< @a; $i++) {
			$POS2QUAL{$chr."#".$pos} = $a[$i];
			$pos++;
		}
	}
}


close QUAL;


print STDERR "Let's go\n";

while (<SNP>) {
	my @a = split " ", $_;

	if (defined( $POS2QUAL{$a[0]."#".$a[1]})) {	
		print $a[0], "\t", $a[1], "\t", $POS2QUAL{$a[0]."#".$a[1]}, "\n";
	}
	else {
		die("Should not happen\n");
	}
}


