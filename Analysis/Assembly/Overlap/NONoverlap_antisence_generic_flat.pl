#! /usr/bin/perl
use strict;

my $usage = "$0 gff segments\n";

my $file1 = shift or die $usage;
my $file2 = shift or die $usage;


### Read gene annotation file (gff)
my %F1=();
open FILE1, $file1 or die "Cannot open annotation file\n";

while (<FILE1>) {
	my @a = split "\t";
	for (my $i = $a[3]; $i<= $a[4]; $i++) {
		$F1{$a[0]."#".$a[6]."#".$i} = 1;
	}
}
close FILE1;


### Read coverage segment file and check overlap with annotation
open FILE2, $file2 or die "Cannot open segments file\n";

while ( my $line = <FILE2>) {

	chomp($line);
	my @a = split("\t", $line);
	my $overlaps = 0;
	
	for (my $i = $a[1]; $i <= ($a[1] + $a[2]); $i++) {
		
		if ( exists $F1{$a[0]."#".$a[3]."#".$i} ) {
			$overlaps++;
			last;
		}
	}
	
	if( $overlaps == 0 && ($a[2] >= 40 || $a[7] >= 20) ) {
		print "$line\n";
	}
}

exit;
