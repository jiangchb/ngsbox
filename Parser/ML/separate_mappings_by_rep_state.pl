#! /usr/bin/perl

use strict;

my $usage = "perl separate_mappings_by_rep_state.pl map.list\n";
my $file = shift or die $usage;
open FILE, $file;

my %FIRST = ();
my %SECOND = ();

while (my $line = <FILE>) {
	my @a = split " ", $line;
	my $rep_state = "";
	if ($a[6] == 1) {
		$rep_state = 1;
	}
	else {
		$rep_state = 0;
	}

	if ($a[9] == 1 or $a[9] == 4 or $a[9] == 5 or $a[9] == 6 or $a[9] == 9 or $a[9] == 10 or$a[9] == 11) {
		$FIRST{$a[3]} = $rep_state;
	}
	else {
		$SECOND{$a[3]} = $rep_state;
	}
}

close FILE;

open FILE, $file;
open REPREP, ">".$file.".rep-rep";
open REPUNIQ, ">".$file.".rep-uniq";
open UNIQUNIQ, ">".$file.".uniq-uniq";

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (defined($FIRST{$a[3]}) and defined($SECOND{$a[3]})) {
		if ($FIRST{$a[3]} == 0 and $SECOND{$a[3]} == 0) {
			print REPREP $line;
		}
		if ($FIRST{$a[3]} == 1 and $SECOND{$a[3]} == 0) {
                        print REPUNIQ $line;
                }
		if ($FIRST{$a[3]} == 0 and $SECOND{$a[3]} == 1) {
                        print REPUNIQ $line;
                }
		if ($FIRST{$a[3]} == 1 and $SECOND{$a[3]} == 1) {
                        print UNIQUNIQ $line;
                }
	}
}




