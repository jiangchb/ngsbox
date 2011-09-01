#!/usr/bin/perl

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;

package Indel;

sub new {
	my $self = {
		ecotype       => '',
		chromosome    => '',
		begin         => 0,
		end           => 0,
		seq           => '',
		support       => 0,
		concordance   => 0,
		repetitive    => 0,
		stype         => 'intergenic',
		gene_id       => 'NA',
		gene_pos      => 0,
		cds_pos       => 0,
		codon_pos     => 0,
		new_stop      => 0,
		lost_stop     => 0,
		splicechange  => 0,
		domain_change => {},
	};
	bless $self;
	return $self;
}

sub init
{
	my ($self, $ecotype, $chromosome, $begin, $end, $seq, $support, $concordance, $repetitive) = @_;
	$self->{ecotype}     = $ecotype;
	$self->{chromosome}  = $chromosome;
	$self->{begin}       = $begin;
	$self->{end}         = $end;
	$self->{support}     = $support;
	$self->{concordance} = $concordance;
	$self->{repetitive}  = $repetitive;
	$self->{seq}         = uc($seq);
}

1;