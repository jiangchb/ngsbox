#! /usr/bin/perl
use strict;

my $usage = "$0 quality_variants unsequenced\n";

open FILE, shift or die $usage;
my $unseq = shift or die $usage;


my $border_size_inner = 50;
my $border_size_outer = 50;

my %MM = ();

while (<FILE>) {
	my @a = split " ";
	$MM{$a[1]."#".$a[2]} = 1;
}

open UNSEQ, $unseq or die $usage;
while(my $line = <UNSEQ>){
	my @a = split " ", $line;
	my $chr = $a[1];
	my $start = $a[2];
	my $end = $a[2];
	my $c_up = 0;
	my $c_down = 0;

	for (my $i = $start + $border_size_inner; $i >= $start - $border_size_outer; $i--) {
		if (defined($MM{$chr."#".$i})) {
			$c_up++;
		}
	}

	for (my $i = $end - $border_size_inner; $i <= $end + $border_size_outer; $i++) {
                if (defined($MM{$chr."#".$i})) {
                        $c_down++;
                }
        }

	#if ($c_up > 0 or $c_down > 0) {
		chomp($line);
		print $line, "\t$c_up\t$c_down\n";
	#}

}




