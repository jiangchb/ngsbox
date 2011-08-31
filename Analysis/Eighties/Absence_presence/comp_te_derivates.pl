#! /usr/bin/perl
use strict;

my $usage = "$0 derivatefile deletionfile percent_overlap\n";

my $file = shift or die $usage;
my $del = shift or die $usage;
my $perc = shift or die $usage;

open FILE, $file;
open DEL, $del;
open CENT, "/ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.centromeres.txt";


########################################################################
# Read in centromeres to see which TEs cannot be assessed for
# absence/presence polymorphism.
########################################################################
my %CEN = ();
while (my $line = <FILE>) {
	my @a = split " ", $line;
	for (my $i = $a[1]; $i <= $a[2]; $i++) {
		$CEN{$a[0]."#".$i} = 1;
	}	
}

########################################################################
# Read in derivatefiles
########################################################################
my %DERIVATE = ();
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (substr($a[0], 0, 1) ne "#") {
		$DERIVATE{$a[0]} = $a[1]."#".$a[2]."#".$a[3];
		for (my $i = $a[2]; $i <= $a[3]; $i++) {
			if (defined($CEN{$a[1]."#".$i})) {
				$i = $a[3]+1; # to break loop
				print STDERR "Excluded from SV analysis:\t", $a[0], "\n";
			}
		}
	}
}

########################################################################
# Read in deletions
########################################################################
my %DELETIONS = ();
my $deleco = "";
while (my $line = <DEL>) {
	my @a = split " ", $line;
	if ($a[0] ne $deleco) {
		$deleco = $a[0];
		#print STDERR "  reading ", $deleco, "\n";
	}

	if (!(defined($DELETIONS{$a[0]}))) {
		%{$DELETIONS{$a[0]}} = ();
	}

	for (my $i = $a[3]; $i <= $a[4]; $i++) {
		${$DELETIONS{$a[0]}}{$a[2]."#".$i} = 1;
		#print STDERR $a[0], "\t", $a[2], "#" , $i, "\n";
	}

	#$DELETIONS{$a[0]} = \%del;
}

########################################################################
# Overlap
########################################################################
open OUT, ">$file.$perc%_deleted"; 
print OUT "#1=present\t0=absent\n";
print OUT "#Accession";
foreach my $derivate (sort keys %DERIVATE) {
	print OUT "\t", $derivate;
}
print OUT "\n";

foreach my $eco (keys %DELETIONS) {
	#print STDERR $eco, "\n";
	my %deletions = %{$DELETIONS{$eco}};
	#foreach my $key (keys %deletions) { print $key, "\t", $deletions{$key}, "\n"; }
	my %deleted_derivates = ();
	print OUT $eco;
	foreach my $derivate (sort keys %DERIVATE) {
		my @a = split "#", $DERIVATE{$derivate};
		my $devri_chr = $a[0];
		my $devri_start = $a[1];
		my $devri_end = $a[2];

#print STDERR "devri: ", $devri_chr, "\t", $devri_start, "\t",  $devri_end, "\n";

		my $del_pos = 0;

		for (my $i = $devri_start; $i <= $devri_end; $i++) {
			if (defined($deletions{$devri_chr."#".$i})) {
				$del_pos++;	
			}
		}

		if (100*($del_pos/($devri_end-$devri_start+1)) >= $perc) {
			print OUT "\t0";
		}
		else {
			print OUT "\t1";
		}
	}
	print OUT "\n";
}


