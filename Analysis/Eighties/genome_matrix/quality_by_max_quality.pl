#! /usr/bin/perl

use strict;

my $usage = "$0 genomeMatrixMaxQual";

my $chr = 1;
my $begin = 1000000;
my $end =   1500000;

my $file = shift or die $usage;

open FILE, $file or die $usage;
open AFILE, "> quality_a.txt" or die $usage;
open CFILE, "> quality_c.txt" or die $usage;
open GFILE, "> quality_g.txt" or die $usage;
open TFILE, "> quality_t.txt" or die $usage;
open DFILE, "> quality_d.txt" or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;

	if ($a[1] % 100000 == 0) {
		print $a[1], "\n";
	}

	if ($a[0] == $chr and $a[1] >= $begin) {

		my $max_a = $a[164];
		my $max_c = $a[166];
		my $max_g = $a[168];
		my $max_t = $a[170];

		for (my $i = 3; $i <= 162; $i+=2) {
			if ($a[$i] eq "A") {
				print AFILE $a[$i+1], "\t", $max_a, "\n";
			}
			elsif ($a[$i] eq "C") {
				print CFILE $a[$i+1], "\t", $max_c, "\n";
			}
			elsif ($a[$i] eq "T") {
				print TFILE $a[$i+1], "\t", $max_t, "\n";
        	        }
			elsif ($a[$i] eq "G") {
				print GFILE $a[$i+1], "\t", $max_g, "\n";
                	}
			elsif ($a[$i] eq "-") {
                                print DFILE $a[$i+1], "\t", $max_g, "\n";
                        }
		}
	}

	if ($a[0] == $chr and $a[1] >= $end) {
		exit(0);
	}

}


