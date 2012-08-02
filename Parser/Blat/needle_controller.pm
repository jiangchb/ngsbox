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
#  Module: Parser::Blat::needle_controller.pm
#  Purpose:
#  In:
#  Out:
#


use FindBin;
use lib $FindBin::Bin;
use needle_alignment;

package needle_controller;


### Constructor
sub new {
	my $self = {
		reads => {},
		reference  => {},
	};
	bless $self;
	return $self;
}


sub align {
	my ($self, $chr, $chrstart, $chrend, $read, $dir) = @_;

	my $align = needleman_wunsch(substr($self->{reference}{$chr}, $chrstart, $chrend-$chrstart+1), $self->{reads}{$read}, $dir);
	#my ($snp_ref, $del_ref, $ins_ref) = parse_alignment($align, 0);

}
	

sub init {
	my ($self, $ref_file, $read_file) = @_;

	$self->parse_fasta(\$self->{reference}, $ref_file);
	$self->parse_fasta(\$self->{reads}, $read_file);

}

sub parse_fasta {

	my ($self, $cont_ref, $file) = @_;

	my $id = "";
	my $seq = "";

	open FASTA, $file or die "Cannot open file ".$file."\n";
	while (my $line = <FASTA>) {
		chomp($line);
		if (substr($line, 0, 1) eq ">") {
			if ($seq ne "") {
				${$cont_ref}{$id} = $seq;
			}
			$seq = "";
			my @a = split " ", $line;
			$id = substr($a[0], 1, length($a[0])-1);
		}
		else {
			$seq .= $line;
		}
	}
	if ($seq ne "") {
		${$cont_ref}{$id} = $seq;
        }

}



1;
