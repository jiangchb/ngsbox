#!/usr/bin/perl

use strict;
use warnings;

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
