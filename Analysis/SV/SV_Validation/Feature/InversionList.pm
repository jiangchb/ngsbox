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
#  Module: Analysis::SV::SV_Validation::Feature::InversionList.pm
#  Purpose:
#  In:
#  Out:
#


use Inversion;

package SNPlist;

sub new {
	my $self = {
		ecotypes   => {},
	};
	bless $self;
	return $self;
}
		

sub get {
	my ($self, $dbh, $chromosome, @ecotypes) = @_;

	my $supress = '';
	if($upper_limit == 0) { $supress = "#" };
	
	### Initialize ecotype SNP hashes
	my %inversions = ();
	$self->{ecotypes}{NR} = \%inversions;
	foreach my $ecotype (@ecotypes) {
		my %inversions = ();
		$self->{ecotypes}{$ecotype} = \%inversions;
	}

	### Select all inversions for ecotypes
	my $q = "";
	
	### Standard quality
	$q = "	SELECT	*
		FROM    sv_inserion
		WHERE	chromosome = $chromosome
		ORDER BY sample, type, chromosome, cluster1_region1_begin";	

	my $sth = $dbh->prepare($q);
	$sth->execute();
	
	### Loop through results
	while(my $ref = $sth->fetchrow_hashref()) {
		my $type        = $ref->{'type'};
		my $ecotype     = $ref->{'sample'};
		my $c1r1_begin  = $ref->{'cluster1_region1_begin'};
		my $c1r1_end    = $ref->{'cluster1_region1_end'};
		my $c1r2_begin  = $ref->{'cluster1_region2_begin'};
		my $c1r2_end    = $ref->{'cluster1_region2_end'};
		my $c2r1_begin  = $ref->{'cluster2_region1_begin'};
                my $c2r1_end    = $ref->{'cluster2_region1_end'};
                my $c2r2_begin  = $ref->{'cluster2_region2_begin'};
                my $c2r2_end    = $ref->{'cluster2_region2_end'};
		my $ori1        = $ref->{'cluster1_orientation'};
		my $ori2        = $ref->{'cluster2_orientation'};
		
		$self->{ecotypes}{$ecotype}{$c1r1_begin} = new Inversion();
		$self->{ecotypes}{$ecotype}{$c1r1_begin}->init($chromosome, $position, $reference, $replacement, $ecotype );
		
		#if( ! exists $self->{ecotypes}{NR}{$c1r1_begin."#".$c2r2_end} ) {
		#	$self->{ecotypes}{NR}{$position} = new Inversion();
		#	$self->{ecotypes}{NR}{$position}->init( $chromosome, $position, $reference, $replacement, $ecotype );
		#}
	}

	return(1);
}

1;
