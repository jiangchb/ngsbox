#!/usr/bin/perl 

use strict;
use warnings;
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
