#!/usr/bin/perl

# --------------------------------------------------------------------
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;

package SNP;

sub new {
	my $self = {
		ecotype       => '',
		chromosome    => 0,
		position      => 0,
		stype         => 'intergenic',
		gene_id       => 'NA',
		gene_pos      => 0,
		cds_pos       => 0,
		codon_pos     => 0,
		ns_change     => 0,
		new_stop      => 0,
		lost_stop     => 0,
		splicechange  => 0,
		ref_base      => '',
		new_base      => '',
		ref_aa        => '',
		new_aa        => '',
		domain_change => {},
	};
	bless $self;
	return $self;
}

sub init
{
	my ($self, $ecotype, $chromosome, $position, $ref_base, $new_base) = @_;
	$self->{ecotype}     = $ecotype;
	$self->{chromosome}  = $chromosome;
	$self->{position}    = $position;
	$self->{ref_base}    = uc($ref_base);
	$self->{new_base}    = uc($new_base);
}

1;
