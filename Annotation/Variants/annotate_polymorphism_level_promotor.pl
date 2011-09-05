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
#  Module: Annotation::Variants::annotate_polymorphism_level_promotor.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 AnnotationGFF quality_reference.txt quality_variants.txt\n";
my $anno = shift or die $usage;
my $reference = shift or die $usage;
my $variants = shift or die $usage;

my %REF = ();
my %POLY = ();
my %PROMOTOR = ();
my %ANNOTATION = ();

###############################################
## Read in files
get_anno_file(\%PROMOTOR, \%ANNOTATION, $anno);
print STDERR "Got anno file\n";
get_qual_file(\%REF, $reference);
print STDERR "Got reference file\n";
get_qual_file(\%POLY, $variants);
print STDERR "Got variants file\n";


GENE: foreach my $tair (sort keys %PROMOTOR) {

	#print STDERR ">$tair<",$GENES{$tair},"\n";

	my @a = split "#", $PROMOTOR{$tair};
	my $chr = $a[0];
	my $begin = $a[1];
	my $end = $a[2];

	###########################################
	# Check for ref call
	my $c = 0;
	my $c_cons = 0;
	my $c_poly = 0;

	for (my $i = $begin; $i <= $end; $i++) {
		$c++;
		# Conservered of polymorph
		if (defined($REF{$chr."#".$i})) {
			$c_cons++;
		}
		elsif (defined($POLY{$chr."#".$i})) {
                        $c_poly++;
                }
	}

	print $tair, "\t", $ANNOTATION{$tair}, "\t", $c, "\t", $c_cons, "\t", $c_poly, "\t";
	if (($c_cons+$c_poly) == 0) {
		print "-1";
	}
	else {
		printf("%.3f", $c_cons/($c_cons+$c_poly));
	}
	print "\n";

}






###################################################
# subroutines

sub get_qual_file {
	my ($ref, $file) = @_;
	
	open FILE, $file or die "Cannot open file: $file\n";
	while (my $line=<FILE>) {
		my @a = split " ", $line;
		if ($a[5] >= 25) {
			${$ref}{$a[1]."#".$a[2]} = $a[4]."#".$a[8];
		}
	}
	close FILE;
}

sub get_anno_file {
	my ($promotor_ref, $anno_ref, $file) = @_;

	open FILE, $file or die "Cannot open file\n";

	while (my $line=<FILE>) {
		my @a = split " ", $line;
		$a[0] =~ s/Chr//g;
		if ($a[2] eq "gene" or $a[2] eq "pseudogene") {
			my @b = split ";", $a[8];
			$b[0] =~ s/ID=//g;
			if ($a[6] eq "+") {
				${$promotor_ref}{$b[0]} = $a[0]."#".($a[3]-1000)."#".($a[3]-1);
			}
			else {
				${$promotor_ref}{$b[0]} = $a[0]."#".($a[4]+1)."#".($a[4]+1000);
			}
			${$anno_ref}{$b[0]} = $a[2];
		}
	}

	close FILE or die "file won't close\n";
}












