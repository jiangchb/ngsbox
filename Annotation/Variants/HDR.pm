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
#  Module: Annotation::Variants::HDR.pm
#  Purpose:
#  In:
#  Out:
#


# --------------------------------------------------------------------
# NGSBox: Variant Annotation
# 
# Annotate SNPs, indels, SVs and CNVs against any eukaryotic gene
# annotation in GFF format
# 
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------


package HDR;

sub new {
	my $self = {
		ecotype         => '',
		chromosome      => '',
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
	my ($self, $ecotype, $chromosome, $begin, $end) = @_;
	$self->{ecotype}     = $ecotype;
	$self->{chromosome}  = $chromosome;
	$self->{begin}       = $begin;
	$self->{end}         = $end;
}

1;
