#!/usr/bin/perl

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;

package HDR;

sub new {
	my $self = {
		ecotype         => '',
		chromosome      => 0,
		begin           => 0,
		end             => 0,
		seq             => '',
		stype           => 'intergenic',
		dtype           => '',
		gene_absence    => {},
		cds_absence     => {},
		utr_absence     => {},
		partial_absence => {}, 
		gene_pos        => 0,
		cds_pos         => 0,
		codon_pos       => 0,
		new_stop        => 0,
		lost_stop       => 0,
		splicechange    => 0,
		domain_change   => {},
	};
	bless $self;
	return $self;
}

sub init
{
	my ($self, $ecotype, $chromosome, $begin, $end, $dtype) = @_;
	$self->{ecotype}     = $ecotype;
	$self->{chromosome}  = $chromosome;
	$self->{begin}       = $begin;
	$self->{end}         = $end;
	$self->{dtype}       = $dtype;
}

1;
