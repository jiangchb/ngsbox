#!/usr/bin/perl

our $VERSION = '1.0';

use strict;
use warnings;

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
