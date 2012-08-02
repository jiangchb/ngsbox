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
#  Module: Parser::Blat::Alignment.pm
#  Purpose:
#  In:
#  Out:
#


package Alignment;


### Constructor
sub new {
	my $self = {
		match  => 0,
		blocks  => 0,
		read_length  => 0,
		read_starts => [],
		read_ends => [],
		read_seqs => [],
		target_id => "",
		target_starts => [],
                target_ends => [],
                target_seqs => [],
	};
	bless $self;
	return $self;
}


sub get_maplist_format {
	my ($self, $hits) = @_;

	return $self->{target_id}."\t";
}


sub init {

	my ($self, $match, $read_length, $read_starts_ref, $read_ends_ref, $read_block_seq_ref, $target_id, $target_starts_ref, $target_ends_ref, $target_block_seq_ref) = @_;
	
	$self->{match} = $match;
	$self->{read_length} = $read_length;
	$self->{read_starts} = @{$read_starts_ref};
	$self->{read_end} = @{$read_ends_ref};
	$self->{read_block_seq} = @{$read_block_seq_ref};
	$self->{target_id} = @{$target_id};
	$self->{target_starts} = @{$target_starts_ref};
        $self->{target_end} = @{$target_ends_ref};
        $self->{target_block_seq} = @{$target_block_seq_ref};	

}


1;
