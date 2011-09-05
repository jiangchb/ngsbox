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
#  Module: Simulation::Simulate_NGS::contamination.pl
#  Purpose:
#  In:
#  Out:
#


our $VERSION = '1.0';


my $contamination_genome = shift;
my $read_length = shift;
my $total_reads = shift;

### Read sequence from E. coli
open CONT, $contamination_genome or die "Cannot open input file $contamination_genome'\n";
my $contamination_seq = '';
while( <CONT> ) {
	chomp;
	if($_ !~ />/) {
		$_ =~ s/N//g;
		$contamination_seq .= $_;
	}
}


### Creates random sequences (contamination) from E.coli
for( my $i = 0; $i < $total_reads; $i++ ) {
	my $pos = int(rand(length($contamination_seq)));
	my $read = substr($contamination_seq, $pos, $read_length);
	print "$pos\t$read\n";
}

exit(0);
