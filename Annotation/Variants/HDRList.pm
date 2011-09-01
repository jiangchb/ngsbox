#!/usr/bin/perl 

# --------------------------------------------------------------------
# NGSBox: Functional Variant Annotation
# 
# Annotate SNPs, indels, SVs and CNVs against any eukaryotic gene
# annotation in GFF format
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
		format     => "",
		hdr_file   => "",
		hdrs       => {},
	};
	bless $self;
	return $self;
}
		

sub get {

	### Init
	my ($self, $hrd_file, $format) = @_;
	
	$self->{format}     = $format;
	$self->{hdr_file}   = $hrd_file;

	open HDRFILE, $hrd_file or die "Cannot open HDR file: $hrd_file\n";


	### SHORE generic HDR format:  sample | chromosome | begin | end | length | ...
	if($format eq "shore") {
		while (<HDRFILE>) {
			chomp;

			my @e = split "\t";

			$self->{hdrs}{$e[1]}{$e[2]} = new HDR;
			$self->{hdrs}{$e[1]}{$e[2]}->init( $e[0], $e[1], $e[2], $e[3]);
		}
	}

	close HDRFILE;

	return(1);
}

1;
