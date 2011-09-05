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
#  Module: Annotation::Variants::SNP.pm
#  Purpose:
#  In:
#  Out:
#


# --------------------------------------------------------------------
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------


package SNP;

sub new {
	my $self = {
		ecotype       => '',
		chromosome    => '',
		position      => 0,
		quality       => 0,
		support       => 0,
		concordance   => 0,
		repetitive    => 0,
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
	my ($self, $ecotype, $chromosome, $position, $ref_base, $new_base, $quality, $support, $concordance, $repetitive) = @_;
	$self->{ecotype}     = $ecotype;
	$self->{chromosome}  = $chromosome;
	$self->{position}    = $position;
	$self->{ref_base}    = uc($ref_base);
	$self->{new_base}    = uc($new_base);
	$self->{quality}     = $quality;
	$self->{support}     = $support;
	$self->{concordance} = $concordance;
	$self->{repetitive}  = $repetitive;
}

1;
