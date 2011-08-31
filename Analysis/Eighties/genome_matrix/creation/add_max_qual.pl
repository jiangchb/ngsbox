#! /usr/bin/perl

use strict;

my $usage = "$0 genomeMatrix\n";

my $file = shift or die $usage;

open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;

	my $max_a = -1;
	my $max_c = -1;
	my $max_g = -1;
	my $max_t = -1;
	my $max_d = -1;

	my $freq_a = 0;
	my $freq_c = 0;
	my $freq_g = 0;
	my $freq_t = 0;
	my $freq_d = 0;
	my $freq_n = 0;

	for (my $i = 3; $i < @a; $i+=2) {
		if ($a[$i] eq "A") {
			if ($a[$i+1] > $max_a) {
				$max_a = $a[$i+1];
			}
			$freq_a++;
		}
		elsif ($a[$i] eq "C") {
			if ($a[$i+1] > $max_c) {
                                $max_c = $a[$i+1];
                        }
                        $freq_c++;
		}
		elsif ($a[$i] eq "T") {
			if ($a[$i+1] > $max_t) {
                                $max_t = $a[$i+1];
                        }
                        $freq_t++;
                }
		elsif ($a[$i] eq "G") {
			if ($a[$i+1] > $max_g) {
                                $max_g = $a[$i+1];
                        }
                        $freq_g++;
                }
		elsif ($a[$i] eq "-") {
                        if ($a[$i+1] > $max_d) {
                                $max_d = $a[$i+1];
                        }
                        $freq_d++;
                }
		elsif ($a[$i] eq "N") {
                        $freq_n++;
                }
	}

	chomp($line);
	print $line, "\t", $freq_a, "\t", $max_a, "\t", $freq_c, "\t", $max_c, "\t", $freq_g, "\t", $max_g, "\t", $freq_t, "\t", $max_t, "\t", $freq_d, "\t", $max_d, "\n";

}


