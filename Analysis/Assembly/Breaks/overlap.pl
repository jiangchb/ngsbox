#! /usr/bin/perl
use strict;
my $usage = "$0 repeatfile locationfile\n";
open REPEAT, shift or die $usage;
open LOC, shift or die $usage;

my %R = ();

while (<REPEAT>) {
	my @a = split " ";
	$R{$a[0]."#".$a[1]} = 1;
}


my $count = 0;
my $count_repeat = 0;

while (<LOC>) {
        my @a = split " ";
	my $begin; my $end;
	if ($a[0] eq "BEG") {
		$begin = $a[2]-($a[3]-1)-200;
		$end = $a[2]-($a[3]-1);
	}
	elsif ($a[0] eq "END") {
		$begin = $a[2]+($a[3]-1);
		$end = $a[2]+($a[3]-1)+200;
	}
	for (my $i = $begin; $i <= $end; $i++) { 
		if (defined($R{$a[1]."#".$i})) {
			$count_repeat++;
		}
		$count++;
	}
}

print "Querried: ", $count, "\n";
print "Repeat: ", $count_repeat, "\n";
print "Perc repeat: ", $count_repeat/$count, "\n";
print "All:  105236502\n";
print "All repeat: 8619290\n";
print "All Perc repeat: ", 8619290/105236502, "\n";
