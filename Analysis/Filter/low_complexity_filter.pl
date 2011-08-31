#!/usr/bin/perl

use strict;
use warnings;

my $file = shift;

open IN, $file or die "Cannot open input file\n";
while( <IN> ) {
	chomp;
	my @elements = split(/\t/, $_);
	if( &seq_complexity($elements[1]) ) {
		print "$_\n";
	}
}
close IN;

exit(0);


sub seq_complexity
{
        my $seq = shift;
	my $length = length($seq);
        my %occ   = (A => 0, C => 0,G => 0,T => 0);
	my %count = (A => 0, C => 0,G => 0,T => 0);

        for(my $i = 0; $i < length($seq); $i++) {
                $count{substr($seq, $i, 1)}++;
		$occ{substr($seq, $i, 1)} = 1;
        }

        if( 
		( ($occ{A} + $occ{C} + $occ{G} + $occ{T}) >= 3) &&
		( ($length - $count{A}) > 4 ) &&
		( ($length - $count{C}) > 4 ) &&
		( ($length - $count{G}) > 4 ) &&
		( ($length - $count{T}) > 4 )
	) {
                return(1);
        }
        else {
                return(0);
        }
}
