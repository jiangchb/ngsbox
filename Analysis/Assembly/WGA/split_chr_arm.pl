#! /usr/bin/perl
use strict;
use warnings;

my $usage = "$0 fastafile\n";

my $file = shift or die $usage;
open FILE, $file or die $usage;

my $seq  = "XX";
my $id   = "";
my $desc = "";

### Parse assembly file
while (my $line = <FILE>) {
	chomp($line);

	# Header found
	if (substr($line, 0, 1) eq ">") {

		# Split header
		my @e = split(" ", $line);

		if ($seq ne "XX") {
			print OUT "$id | $desc\n$seq\n";

			# Open new file if new chromosome arm is found
			if($e[2] ne $desc) {
				close OUT;
				$desc = $e[2];
				open OUT, ">$desc.fa" or die "Cannot open $file.$desc.fa\n";
			}
		}
		else {
			$desc = $e[2];
			open OUT, ">$desc.fa" or die "Cannot open $file.$desc.fa\n";
		}

		# Reset
		$seq  = "";
		$id   = $e[0];
		$desc = $e[2];
	}
	else {
		$seq .= $line;
	}
}

### Last one
if ($seq ne "") {
	print OUT "$id | $desc\n$seq\n";
}

	

