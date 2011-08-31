#!/usr/bin/perl 

use strict;
use warnings;

package Mutations;

sub new {
	my ($self) = @_;
	$self = {
		chromosome	=> (),
		position	=> (),
		sample	    	=> (),
		fraction     	=> (),
	};
	bless $self;
	return $self;
}

#############################################################################################
sub set {
	my ($self) = @_;	

	$self->{chromosome}[0] = 4;
	$self->{position}[0] = 16702262;
	$self->{sample}[0] = "Col-0";
	$self->{fraction}[0] = 1;

	return(1);
}


1;
