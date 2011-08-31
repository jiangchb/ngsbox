#! /usr/bin/perl
use strict;

my %R = ();

while (<STDIN>) {
	my ($chr, $pos, $alg, $id) = split " ";	
	$R{$id} .= $chr."#";
}

my $c = 0;
my $c_g = 0;
my $c_o = 0;

foreach my $val (values %R) {
	my @chrs = split "#", $val;
	my $o_flag = 0;
	my $g_flag = 0;
	foreach my $i (@chrs) {
		if ($i==702 || $i==696 || $i==697 || $i==698 || $i==699 || $i==700 || $i==701) {
			$o_flag = 1;
		} else {
			$g_flag = 1;
		}
	}	
	# set counts
	if ($o_flag == 0) {
		$c_g++;
	} elsif ($g_flag == 0) {
		$c_o++;
	}
	$c++;
}

print STDOUT "all:\t\t", $c, "\n";
print STDOUT "ambi:\t\t", $c-($c_g+$c_o), "\n";
print STDOUT "genomic:\t", $c_g, "\n";
print STDOUT "organelles:\t", $c_o, "\n";



