#! /usr/bin/perl

use strict;

# Makes a "reference sequence" with inserted snps.

my $usage = "make_ref_with_inserted_snps.pl snpfile fastafile\n";
my $snpfile = shift or die $usage;
my $fastafile = shift or die $usage;

open SNP, $snpfile;

my %SNPS = ();

while(my $line = <SNP>) {
	my @a = split " ", $line;
	$SNPS{$a[1]."#".$a[2]} = $a[4];
}

close SNPS;

open FASTA, $fastafile;

my $chr = 0;
my $pos = 0;
while (my $line = <FASTA>) {
	chomp($line);
	if (substr($line, 0, 1) eq ">") {
		$chr++;
		$pos = 0;
		print $line, "\n";
	}
	else {
		for (my $i = 0; $i < length($line); $i++) {
			$pos++;
			if (defined($SNPS{$chr."#".$pos})) {
				print $SNPS{$chr."#".$pos};
			}
			else {
				print substr($line, $i, 1);
			}
		}
		print "\n";
	}
}

