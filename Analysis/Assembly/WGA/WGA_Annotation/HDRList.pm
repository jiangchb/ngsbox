#!/usr/bin/perl 

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;
use HDR;

package HDRList;

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

	### Shore assembly HDR format:  sample | type | chromosome | scaffold | HSP number | r-begin | r-end | s-begin | s-end | r-length | s-length | seq
	while (<INDELFILE>) {
		chomp;
		
		my @e = split "\t";

		if ( $e[2] == $chromosome ) {
			$self->{indels}{$e[5]} = new HDR;
			$self->{indels}{$e[5]}->init( $e[0], $e[2], $e[5], $e[6], $e[1]);
		}
	}

	close INDELFILE;

	return(1);
}

1;
