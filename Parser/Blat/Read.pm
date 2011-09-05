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
#  Module: Parser::Blat::Read.pm
#  Purpose:
#  In:
#  Out:
#


use Alignment;

package Read;

### Constructor
sub new {
	my $self = {
		best_match  => 0,
		type  => "",
		alignments  => [],
	};
	bless $self;
	return $self;
}
		
sub set_deletion_alignment {
	my ($self) = @_;

	for (my $i = 0; $i < @{$self->{alignments}}; $i++) {
		
        }	

}


sub get_maplist_format {
	my ($self) = @_;

	my $maplist = "";
	
	foreach my $align (@{$self->{alignments}}) {
		$maplist .= $align->get_maplist_format(@{$self->{alignments}}+0);
	}

}

sub add_alignment {

	my ($self, $match, $c_read_length, $c_read_starts_ref, $c_read_ends_ref, $c_read_block_seq_ref, $c_target_id, $c_target_starts_ref, $c_target_ends_ref, $c_target_block_seq_ref) = @_;

	if ($self->{best_match} < $match) {
		$self->{best_match} = $match;
	}

	my $alignment = new Alignment();
	$alignment->init($match, $c_read_length, $c_read_starts_ref, $c_read_ends_ref, $c_read_block_seq_ref, $c_target_id, $c_target_starts_ref, $c_target_ends_ref, $c_target_block_seq_ref);

	push @{$self->{alignments}}, $alignment;

}

sub remove_repetitive_alignments {
	my ($self, $weak) = @_;

	for (my $i = 0; $i < @{$self->{alignments}}; $i++) {
		if (${$self->{alignments}}[$i]->{match} < $self->{best_match} * (1-($weak/100))) {
			delete ${$self->{alignments}}[$i];
			$i--; 
		}
	}

}


1;
