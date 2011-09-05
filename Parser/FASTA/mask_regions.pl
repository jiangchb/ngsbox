#! /usr/bin/perl
use strict;
use warnings;

###### 
# NGSbox - bioinformatics analysis tools for next generation sequencing data
#
# Copyright 2007-2011 Stephan Ossowski, Korbinian Schneeberger
# 
# NGSbox is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or any later version.
#
# NGSbox is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# Please find the GNU General Public License at <http://www.gnu.org/licenses/>.
#
#  -------------------------------------------------------------------------
#
#  Module: Parser::FASTA::mask_regions.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 file chr start end\n";

my $file = shift or die $usage;
my $chr = shift or die $usage;
my $pos = shift or die $usage;
my $end = shift or die $usage;

open FILE, $file or die $usage;

my $flag = 0;
my $id = "";
my $seq = "";
while (my $line = <FILE>) {
	chomp($line);
	if (substr($line, 0, 1) eq ">") {
		print_mask($seq, $id, $chr, $pos, $end);
		my @a = split " ", $line;
		$id = substr($a[0], 1, length($a[0])-1);
		$seq = "";
	}
	else {
		$seq .= $line;
	}	
}
print_mask($seq, $id, $chr, $pos, $end);


sub print_mask {
	my ($seq, $id, $chr, $pos, $end) = @_;


	if ($id ne "") {
		print ">", $id, "\n";
		my $line_break = 0;
	        if ($chr == $id) {
			my $lb = 0;
			for (my $i = 0; $i<$pos-1; $i++) {
	        		print substr($seq, $i, 1);
				$lb++;
				if ($lb%79 == 0) {
					print "\n";
				}
			}
                	for (my $i = 0; $i < $end-$pos+1; $i++) {
                		print "N";
				$lb++;
				if ($lb%79 == 0) {
                                        print "\n";
                                }
	                }
			for (my $i = $end; $i<length($seq); $i++) {
	        	      	print substr($seq, $i, 1);
				$lb++;
				if ($lb%79 == 0) {
                                        print "\n";
                                }
                        }
			if ($lb%79 != 0) {
                        	print "\n";
                        }
	        }
        	else {
			my $i;
			for ($i = 0; $i<length($seq); $i++) {
        			print substr($seq, $i, 1); 
				if (($i+1)%79 == 0) {
					print "\n";
				}
			}
			if ($i+1%79 != 0) {
                                print "\n";
                        }
	        }
	}
}




