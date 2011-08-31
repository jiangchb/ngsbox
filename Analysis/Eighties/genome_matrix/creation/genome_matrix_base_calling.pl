#! /usr/bin/perl

use strict;

my $usage = "$0 genomeMatrixMaxQual MinMaxQual MinQual\n";

my $file = shift or die $usage;
my $MinMaxQual = shift or die $usage;
my $MinQual = shift or die $usage;

open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;

	my $max_a = $a[164];
	my $max_c = $a[166];
	my $max_g = $a[168];
	my $max_t = $a[170];
	my $max_d = $a[172];

#print STDERR $max_a, "\t", $max_c, "\t", $max_g, "\t", $max_t, "\n";

	print $a[0], "\t", $a[1], "\t", $a[2];

	for (my $i = 3; $i <= 162; $i+=2) {
		if ($a[$i] eq "A") {
			if (($max_a >= $MinMaxQual and $a[$i+1] > $MinQual) or ($a[$i+1] > $MinMaxQual)) {
				print "\tA";
			}
			else {
				print "\tN";
			}
		}
		elsif ($a[$i] eq "C") {
			if (($max_c >= $MinMaxQual and $a[$i+1] > $MinQual) or ($a[$i+1] > $MinMaxQual)) {
	                        print "\tC";
			}
                	else {
	                        print "\tN";
			}
                }
		elsif ($a[$i] eq "G") {
			if (($max_g >= $MinMaxQual and $a[$i+1] > $MinQual) or ($a[$i+1] > $MinMaxQual)) {
	                        print "\tG";
	                }
        	        else {
                	        print "\tN";
			}
                }
		elsif ($a[$i] eq "T") {
			if (($max_t >= $MinMaxQual and $a[$i+1] > $MinQual) or ($a[$i+1] > $MinMaxQual)) {
	                        print "\tT";
	                }
        	        else {
                	        print "\tN";
			}
                }
		elsif ($a[$i] eq "-") {
                        if (($max_d >= $MinMaxQual and $a[$i+1] > $MinQual) or ($a[$i+1] > $MinMaxQual)) {
                                print "\t-";
                        }
                        else {
                                print "\tN";
                        }
                }
		#if ($a[$i] ne "A" and $a[$i] ne "C" and $a[$i] ne "G" and $a[$i] ne "T" and $a[$i] ne "-") {
		else {
			print "\tN";
		}
	}
	
	
	#print "\t", $a[163], "\t", $a[164], "\t", $a[165], "\t", $a[166], "\t", $a[167], "\t", $a[168], "\t", $a[169], "\t", $a[170], "\n";
	print "\n";

}


