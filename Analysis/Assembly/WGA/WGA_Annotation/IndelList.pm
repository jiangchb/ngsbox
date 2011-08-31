#!/usr/bin/perl 

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;
use Indel;

package IndelList;

### Constructor
sub new {
	my $self = {
		indel_file   => "",
		chromosome   => "",
		indels       => {},
	};
	bless $self;
	return $self;
}
		

sub get {

	### Init
	my ($self, $chromosome, $indel_file) = @_;
	$self->{chromosome}   = $chromosome;
	$self->{indel_file}   = $indel_file;

	open INDELFILE, $indel_file or die "Cannot open indel file: " . $indel_file . "\n";

	### Shore indel format:  sample | chromosome | begin | end | scaffold | s-begin | s-end | length | seq
	while (<INDELFILE>) {
		chomp;
		
		my ($sample, $chr, $begin, $end, $scaff, $s_begin, $s_end, $length, $seq) = split "\t";

		if ( $chr == $chromosome ) {
			$self->{indels}{$begin} = new Indel;
			$self->{indels}{$begin}->init( $sample, $chromosome, $begin, $end, $seq);
		}
	}

	close INDELFILE;

	return(1);
}

1;
