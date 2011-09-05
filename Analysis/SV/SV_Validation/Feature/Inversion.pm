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
#  Module: Analysis::SV::SV_Validation::Feature::Inversion.pm
#  Purpose:
#  In:
#  Out:
#



package SNP;

sub new {
	my $self = {
		type          => '',
		sample        => '',
		chromosome    => 0,
		cluster1_region1_begin => 0,
		cluster1_region1_end => 0,
		cluster1_region2_begin => 0,
                cluster1_region2_end => 0,
		cluster2_region1_begin => 0,
                cluster2_region1_end => 0,
                cluster2_region2_begin => 0,
                cluster2_region2_end => 0,
		cluster1_orientation => '',
		cluster2_orientation => '',
		gap1          => 0,
		uniq_coverage_gap1 => 0,
		core_uniq_coverage_gap1 => 0,
		rep_coverage_gap1 => 0,
		gap2          => 0,
                uniq_coverage_gap2 => 0,
                core_uniq_coverage_gap2 => 0,
                rep_coverage_gap2 => 0
	};
	bless $self;
	return $self;
}

sub init
{
	my ($self, $type, $sample, $chromosome, ) = @_;
	$self->{type}      = $type;
	$self->{sample}      = $sample;
	$self->{chromosome} = $chromosome;

}

1;
