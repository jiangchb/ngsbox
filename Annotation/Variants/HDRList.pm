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
#  Module: Annotation::Variants::HDRList.pm
#  Purpose:
#  In:
#  Out:
#


# --------------------------------------------------------------------
# NGSBox: Functional Variant Annotation
# 
# Annotate SNPs, indels, SVs and CNVs against any eukaryotic gene
# annotation in GFF format
# 
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use HDR;

package HDRList;

### Constructor
sub new {
	my $self = {
		format     => "",
		hdr_file   => "",
		hdrs       => {},
	};
	bless $self;
	return $self;
}
		

sub get {

	### Init
	my ($self, $hrd_file, $format) = @_;
	
	$self->{format}     = $format;
	$self->{hdr_file}   = $hrd_file;

	open HDRFILE, $hrd_file or die "Cannot open HDR file: $hrd_file\n";


	### SHORE generic HDR format:  sample | chromosome | begin | end | length | ...
	if($format eq "shore") {
		while (<HDRFILE>) {
			chomp;

			my @e = split "\t";

			$self->{hdrs}{$e[1]}{$e[2]} = new HDR;
			$self->{hdrs}{$e[1]}{$e[2]}->init( $e[0], $e[1], $e[2], $e[3]);
		}
	}

	close HDRFILE;

	return(1);
}

1;
