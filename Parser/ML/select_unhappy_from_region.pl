#! /usr/bin/perl


my $usage = "\n$0 map.list chr begin end\n\n";
my $file = shift or die $usage;
my $chr = shift or die $usage;
my $begin = shift or die $usage;
my $end = shift or die $usage;

open FILE, $file or die "Cannot open file\n";

my $written = 0;

while (my $line = <FILE>) {
	my @a = split " ", $line;

	if ($a[0] == $chr && $a[1] >= $begin && $a[1] <= $end) {
		my $pe = ($a[9]-2) % 3;
		#if ($pe == 2) {
		if ($a[9] == 4 || $a[9] == 7 || $a[9] == 10 || $a[9] == 13) {
			print $line;
		}
		$written = 1;
	}
	else {
		if ($written == 1) {
			exit(0);
		}
	}
}


