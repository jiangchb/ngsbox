#! /usr/bin/perl


my $usage = "$0 irfoutput";
my $file = shift or die $usage;

open FILE, $file or die $usage;

my $chr = -1;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (@a > 15) {
		print $chr;
		for (my $i = 0; $i < @a; $i++) {
			print "\t", $a[$i];
		}
		print "\n";
	}
	else {
		if (substr($line, 0, 8) eq "Sequence") {
			$chr = @a[1];
print STDERR "Now parsing chromsome:", $chr, "\n";
		}
	}
}

