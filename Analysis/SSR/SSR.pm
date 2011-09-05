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
#  Module: Analysis::SSR::SSR.pm
#  Purpose:
#  In:
#  Out:
#



package SSR;

sub new {
	my $self = {
		chromosome    => 0,
		position      => 0,
		end           => 0,
		len	      => 0,
		instances     => 0,
		sequence      => '',
		spanning_ungapped      => 0,
		spanning_gapped      => 0,
		variation     => '',
	};
	bless $self;
	return $self;
}

sub set
{
	my ($self, $chromosome, $position, $end, $len, $instance, $sequence) = @_;
	$self->{chromosome} = $chromosome;
	$self->{position}   = $position;
	$self->{end}        = $end;
	$self->{len}        = $len;
	$self->{instances}  = $instance;
	$self->{sequence}   = $sequence;
}

1;
