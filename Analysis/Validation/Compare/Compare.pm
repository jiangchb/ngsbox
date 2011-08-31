#!/usr/bin/perl 

use strict;
use warnings;

package Compare;

sub new {
	my $self = {
	};
	bless $self;
	return $self;
}

#############################################################################################
### Calculates the genome sequence of a sampled genome
sub compare {
	my ($self, $seq1, $seq2) = @_;
	
	my $len_diff = abs(length($seq1) - length($seq2));

	for (my $i = 0; $i <= $len_diff; $i++) {
		my $start1 = 0;
		my $start2 = 0;
		my $len;
		if (length($seq1) > length($seq2)) {
			$start1 = $i;
			$len = length($seq2);
		}
		elsif (length($seq1) < length($seq2)) {
			$start2 = $i;
			$len = length($seq1);
		}
		else {
			$len = length($seq1);
		}
		
		my $ident = 1;

		for (my $j = 0; $j < $len; $j++) {
			if (	substr($seq1, $start1+$j, 1) ne substr($seq2, $start2+$j, 1) and 
				substr($seq1, $start1+$j, 1) ne "N" and
				substr($seq2, $start2+$j, 1) ne "N"
			) {
				$ident = 0;
				$j = $len;
			}
		}

		if ($ident == 1) {
			return 1;
		}
	}

	return 0;
}

1;

