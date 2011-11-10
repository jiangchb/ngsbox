#!/usr/bin/perl

use strict;

open(IN, "<$gencode_file") or die "Can't open $gencode_file.\n";
my %all_genes;
while(<in>){
	next if(/^##/); #ignore header
	chomp;

	my %attribs = ();
	my ($chr, $source, $type, $start, $end, $score, 
		$strand, $phase, $attributes) = split("\t");

	#store nine columns in hash
	my %fields = (
		chr        => $chr,
		source     => $source,
		type       => $type,
		start      => $start,
		end        => $end,
		score      => $score,
		strand     => $strand,
		phase      => $phase,
		attributes => $attributes,
	);

	my @add_attributes = split(";", $attributes);

	# store ids and additional information in second hash
	foreach my $attr ( @add_attributes ) {
		next unless $attr =~ /^\s*(.+)\s(.+)$/;
		$c_type  = $1;
		$c_value = $2;
		if($c_type  && $c_value){
			if(!exists($attribs{$c_type})){
				$attribs{$c_type} = [];
			}
			push(@{ $attribs{$c_type} }, $c_value);
		}
	}

	#work with the information from the two hashes...
	#eg. store them in a hash of arrays by gene_id:
	if(!exists($all_genes{$attribs{'gene_id'}->[0]})){
		$all_genes{$attribs{'gene_id'}->[0]} = [];
	}
	push(@{ $all_genes{$attribs{'gene_id'}->[0]} }, \%fields);
}

print "Example entry ENSG00000223972: ".
	$all_genes->{"ENSG00000223972"}->[0]->{"type"}.", ".
	$all_genes->{"ENSG00000223972"}->[0]->{"chrom"}." ".
	$all_genes->{"ENSG00000223972"}->[0]->{"start"}."-".
	$all_genes->{"ENSG00000223972"}->[0]->{"end"}."\n";
