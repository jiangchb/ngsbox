#! /usr/bin/perl
use strict;
my $usage = "$0 annotationGFF refcalls sample\n";
my $annofile = shift or die $usage;
my $reffile = shift or die $usage;
my $sample = shift or die $usage;

##### Read in ref calls
print STDERR "Start reading refcalls...\n";
my %REF = ();
open REF, $reffile or die "Cannot open file: ", $reffile, "\n";
my $c = 0;
while (my $line = <REF>) {
	$c++;
	if ($c % 1000000 == 0) { print STDERR $line; }
	my @a = split " ", $line;
	if ($a[5] >= 25) {
		$REF{$a[1]."#".$a[2]} = 1;
	}
}
close REF;
print STDERR "...done\n";


open ANNO, $annofile or die "Cannot open file: ", $annofile, "\n";
while (<ANNO>) {
	my @a = split " ";
	if ($a[2] eq "gene") {
		my @b = split ";", $a[8];
		my $id = substr($b[0], 3, length($b[0]) -3);
		my $chr = substr($a[0], 3, 1);
		my $begin; my $end;
		if ($a[6] eq "+") {
			$begin = $a[3] - 2000;
			$end = $a[3] - 1;
		}
		else {
			$begin = $a[4] + 1;
			$end = $a[4] + 2000;
		}
		my $refcount = 0;
		for (my $i = $begin; $i <= $end; $i++) {
			if (defined($REF{$chr."#".$i})) {
				$refcount++;
			} 
		}
		print $sample, "\t", $id, "\t", $refcount, "\n";
	}
}

