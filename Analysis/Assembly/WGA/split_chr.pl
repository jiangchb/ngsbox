#! /usr/bin/perl
use strict;
use warnings;

my $usage = "$0 fastafile\n";

my $file = shift or die $usage;
open FILE, $file or die $usage;

my $seq  = "XX";
my $id   = "";
my $desc = "";
my $chr  = "";

### Parse assembly file
while (my $line = <FILE>) {
	chomp($line);

	# Header found
	if (substr($line, 0, 1) eq ">") {
		
		# Split header
		my @e = split(" ", $line);
		my @f = split("_", $e[2]); 

		if ($seq ne "XX") {
			print OUT "$id | $desc\n$seq\n";

			# Open new file if new chromosome arm is found
			if($f[0] ne $chr) {
				close OUT;
				$chr = $f[0];
				open OUT, ">$chr.fa" or die "Cannot open $chr.fa\n";
			}
		}
		else {
			$chr = $f[0];
			open OUT, ">$chr.fa" or die "Cannot open $chr.fa\n";
		}

		# Reset
		$seq  = "";
		$id   = $e[0];
		$desc = $e[2];
		$chr  = $f[0];
	}
	else {
		$seq .= $line;
	}
}

### Last one
if ($seq ne "") {
	print OUT "$id | $desc\n$seq\n";
}

	

