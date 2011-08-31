#!/usr/bin/perl

use strict;
use warnings;

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
