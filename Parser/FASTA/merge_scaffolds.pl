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
#  Module: Parser::FASTA::fasta_oneliner.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 scaff_file fasta_file\n";

my $scaff = shift or die $usage;
my $fasta  = shift or die $usage;


### Create spacer
my $spacer = "";
for( my $i = 0; $i < 50000; $i++ ) {
	$spacer .= "N";
}

print STDERR length($spacer) . "\n";

### Read fasta file
open FASTA, $fasta or die $usage;

my $seq = "";
my $id  = "";
my %genome = ();
my $counter = 1;

while (my $line = <FASTA>) {
	chomp($line);

	if (substr($line, 0, 1) eq ">") {

		if ($seq ne "") {
			print STDERR "$counter\t$id\n";
			my @tmp = ($counter, $seq);
			$genome{$id} = \@tmp;
			$counter++;
		}
		$seq = "";
		$id = substr($line, 1);
	}
	else {
		$seq .= $line;
	}
}

if ($seq ne "") {
	print STDERR "$counter\t$id\n\n\n";
	my @tmp = ($counter, $seq);
	$genome{$id} = \@tmp;
}


### Read scaffold file
open SCAFF, $scaff or die $usage;

# first scaffold to be merged
my $first = <SCAFF>;
chomp $first;
print STDERR "$first\n\n";
my $merge_counter = $genome{$first}[0];
print STDERR "$merge_counter\n\n";
my $new_id = $first;
print STDERR "$new_id\n\n";
my $newseq = $genome{$first}[1];
delete $genome{$first};

# concatenate the other merge-scaffolds
while($id = <SCAFF>) {
	chomp $id;
	print STDERR "$id\n\n";
	$new_id .= $id;
	print STDERR "$new_id\n\n";
	$newseq = $newseq . $spacer . $genome{$id}[1];
	delete $genome{$id};
}

my @tmp = ($merge_counter, $newseq);
$genome{$new_id} = \@tmp;


foreach my $scaffid (sort {$genome{$a}[0]<=>$genome{$b}[0]} keys %genome) {
	print ">" . $scaffid . "\n" . $genome{$scaffid}[1] . "\n";
}

exit(0);

