#! /usr/bin/perl
use strict;
my $usage = "$0 NBSLRR TAIR8_genes_trasposons shorecout\n";

my $nbslrr = shift or die $usage;
my $tair8 = shift or die $usage;
my $file = shift or die $usage;

#########################################################################################
my %NBS = ();
open FILE, $nbslrr or die "Cannot open file\n";
while (my $line = <FILE>) {
	chomp($line);
	$NBS{$line} = 1;
}
close FILE;

#########################################################################################
my %POS2ANNO = ();
open FILE, $tair8 or die "Cannot open file\n";
while (my $line = <FILE>) {
	if (substr($line, 0, 1) ne "#") {
		my @a = split " ", $line;
		$POS2ANNO{$a[1]."#".$a[2]."#".$a[3]} = $a[4];
	}
}
close FILE;

#########################################################################################
open FILE, $file or die "Cannot open file\n";

while (my $line = <FILE>) {
	if (not(substr($line, 0, 1) eq "#") and length($line) > 5)  {
		my @a = split " ", $line;

		chomp($line);
		print $line;

		if (defined($POS2ANNO{$a[1]."#".$a[2]."#".$a[3]})) {
			print "\t", $POS2ANNO{$a[1]."#".$a[2]."#".$a[3]};
			if (defined($NBS{$POS2ANNO{$a[1]."#".$a[2]."#".$a[3]}})) {
				print "\tNBS_LRR";
			}
			else {
				print "\tNULL";
			}
		}
		else {
			print "\tNULL\tNULL";
		}

		print "\n";
	}
}	

close FILE;



