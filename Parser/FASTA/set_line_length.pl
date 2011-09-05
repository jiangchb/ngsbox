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
#  Module: Parser::FASTA::set_line_length.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 fasta linelength\n";
my $file = shift or die $usage;
my $ll = shift or die $usage;
open FILE, $file or die $usage;

my $seq = "";

while (my $line = <FILE>) {
        chomp($line);
        if (substr($line, 0, 1) eq ">") {
                if ($seq ne "") {
			print_seq();
                }
		print $line, "\n";
		$seq = "";
        }
        else {
                $seq .= $line;
        }
}
print_seq();


sub print_seq {
	my $count = 0;
	for (my $i = 0; $i < length($seq); $i++) {
		$count++;
		print substr($seq, $i, 1);
		print "\n" if $count%$ll == 0;
	}
	print "\n" if $count%$ll != 0;

}


