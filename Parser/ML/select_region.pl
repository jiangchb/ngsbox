#! /usr/bin/perl


my $usage = "\nselect_region.pl map.list chr begin end\n\n";
my $file = shift or die $usage;
my $chr = shift or die $usage;
my $begin = shift or die $usage;
my $end = shift or die $usage;

open FILE, $file or die "Cannot open file\n";
open OUT, ">".$file.".$chr:$begin..$end" or die "Cannot open file\n";

my $written = 0;

while (my $line = <FILE>) {
	my @a = split " ", $line;

	if ($a[1]%200000 == 0) {
		print STDERR $a[0], "\t", $a[1], "\n";
	}

	if ($a[0] == $chr && $a[1] >= $begin && $a[1] <= $end) {
		print OUT $line;
		$written = 1;
	}
	else {
		if ($written == 1) {
			exit(0);
		}
	}
}


